#!/usr/bin/env python3
"""
BioLink — связывающий слой над UniProt, NCBI и RCSB PDB.

Принимает на вход имя гена/белка (например "insulin", "GFP") или UniProt-код
(например "P01308") и собирает данные из трёх независимых открытых баз данных
в одну общую "карточку" объекта:

  1. UniProt      — что за белок, какая у него функция, какой организм.
  2. NCBI (gene)  — данные про сам ген (хромосома, где лежит в геноме, описание).
  3. RCSB PDB     — есть ли расшифрованная 3D-структура этого белка.

Это НЕ система предсказаний — это только "склейка" данных по одному запросу.
Карточка честно показывает, насколько уверенно запрос сопоставлен с записью
(блок "УВЕРЕННОСТЬ"). Подробности — в README.md рядом с этим файлом.

Запуск:
    python3 prototype.py insulin
    python3 prototype.py P01308
    python3 prototype.py GFP
"""

import argparse
import csv
import os
import re
import sys
import time
import json
from datetime import datetime, timezone

import requests

USER_AGENT = "biolink/1.0"
HEADERS = {"User-Agent": USER_AGENT}
TIMEOUT = 20

# Одна сессия на все запросы — переиспользует TCP/TLS-соединение к серверу,
# вместо установки нового соединения на каждый запрос.
SESSION = requests.Session()
SESSION.headers.update(HEADERS)

# Имя файла карточки строится из пользовательского запроса (когда нет
# accession) — ограничиваем длину, чтобы не превысить лимит файловой системы.
MAX_SAFE_NAME_LEN = 60

# У NCBI eutils без API-ключа лимит ~3 запроса/сек — нужна пауза между
# последовательными обращениями.
NCBI_DELAY = 0.4

# NCBI просит идентифицировать инструмент в E-utilities (tool/email).
NCBI_TOOL = "biolink"
NCBI_EMAIL = os.getenv("NCBI_EMAIL")

# Бесплатный API-ключ NCBI (аккаунт NCBI -> Settings -> API Key Management)
# поднимает лимит с ~3 до ~10 запросов/сек. Читается из переменной окружения
# или назначается веб-интерфейсом из .streamlit/secrets.toml.
NCBI_API_KEY = os.getenv("NCBI_API_KEY")

# Пауза между последовательными обращениями к eutils: без ключа — не чаще
# ~3 запросов/сек, с ключом — до 10 (0.15 с оставляет запас).
NCBI_DELAY_WITH_KEY = 0.15


def ncbi_delay():
    """Пауза между запросами к NCBI — зависит от наличия API-ключа."""
    return NCBI_DELAY_WITH_KEY if NCBI_API_KEY else NCBI_DELAY


def ncbi_params(**kwargs):
    params = {"retmode": "json", "tool": NCBI_TOOL}
    if NCBI_EMAIL:
        params["email"] = NCBI_EMAIL
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    params.update(kwargs)
    return params

# UniProt accession имеет два формата: коды на O/P/Q устроены иначе,
# чем на остальных буквах (P01308, Q9Y5Q6, A0A0C5B5G6...).
UNIPROT_ACCESSION_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9]){0,1})$",
    re.IGNORECASE,
)

# Русские бытовые названия распространённых организмов — чтобы в карточке
# писать "Человек (Homo sapiens)": сначала бытовое название, а научное
# (латинское, видовое) — в скобках. Ключи — научные названия из UniProt
# (поле organism.scientificName).
ORGANISM_RU = {
    "Homo sapiens": "Человек",
    "Mus musculus": "Домовая мышь",
    "Rattus norvegicus": "Серая крыса",
    "Cricetulus griseus": "Китайский хомячок",
    "Mesocricetus auratus": "Золотистый хомячок",
    "Oryctolagus cuniculus": "Кролик",
    "Ovis aries": "Овца",
    "Bos taurus": "Домашний бык",
    "Sus scrofa": "Свинья",
    "Gallus gallus": "Курица",
    "Xenopus laevis": "Шпорцевая лягушка",
    "Danio rerio": "Данио-рерио",
    "Drosophila melanogaster": "Плодовая мушка",
    "Anopheles gambiae": "Малярийный комар",
    "Caenorhabditis elegans": "Круглый червь",
    "Saccharomyces cerevisiae": "Пекарские дрожжи",
    "Schizosaccharomyces pombe": "Делящиеся дрожжи",
    "Escherichia coli": "Кишечная палочка",
    "Bacillus subtilis": "Сенная палочка",
    "Pseudomonas aeruginosa": "Синегнойная палочка",
    "Mycobacterium tuberculosis": "Туберкулёзная палочка",
    "Arabidopsis thaliana": "Резуховидка Таля",
    "Oryza sativa": "Рис",
    "Zea mays": "Кукуруза",
    "Nicotiana tabacum": "Табак",
    "Glycine max": "Соя",
    "Aequorea victoria": "Медуза эквория",
    "Macaca mulatta": "Макака-резус",
    "Plasmodium falciparum": "Малярийный плазмодий",
}


def organism_display(scientific_name):
    """Бытовое название организма, научное — в скобках: "Человек (Homo sapiens)".

    Если бытового названия в словаре нет, возвращается только научное
    (латинское) название — оно само по себе является научным названием вида.
    """
    if not scientific_name:
        return "(организм не указан)"
    ru = ORGANISM_RU.get(scientific_name)
    return f"{ru} ({scientific_name})" if ru else scientific_name


# Статус проверки записи UniProt (entryType): прошёл ли он ручную проверку
# кураторов. Reviewed-записи (Swiss-Prot) проверены экспертами, unreviewed
# (TrEMBL) — автоматическая аннотация без ручной проверки.
ENTRY_STATUS_RU = {
    "UniProtKB reviewed (Swiss-Prot)": "Проверено экспертами UniProt (Swiss-Prot)",
    "UniProtKB unreviewed (TrEMBL)": "Не проверено экспертами (TrEMBL)",
}


def entry_status_display(entry_type):
    """Статус проверки записи UniProt по-русски; если статус незнаком — оригинал."""
    return ENTRY_STATUS_RU.get(entry_type, entry_type or "—")


def gene_primary_matches(entry, term):
    # gene_exact матчит не только основное имя гена, но и синонимы
    # (например, у рисового гена NAL1 среди синонимов есть "GFP") —
    # поэтому отдельно проверяем совпадение именно с главным именем.
    for g in entry.get("genes") or []:
        if not g:
            continue
        main_name = (g.get("geneName") or {}).get("value")
        if main_name and main_name.lower() == term.lower():
            return True
    return False


def gene_synonym_matches(entry, term):
    # Общепринятое название может быть только синонимом ("p53" — синоним,
    # главный ген — TP53): тогда это единственный способ найти запись,
    # не откатываясь на неконтролируемый текстовый поиск.
    for g in entry.get("genes") or []:
        if not g:
            continue
        for syn in (g.get("synonyms") or []):
            value = (syn or {}).get("value")
            if value and value.lower() == term.lower():
                return True
    return False


def is_human(entry):
    return (entry.get("organism") or {}).get("taxonId") == 9606


def log(msg):
    print(msg)


def safe_get(url, params=None, retries=2, with_headers=False):
    """GET-запрос с ретраями на временные сбои и разбором JSON.

    Возвращает (data, None), а с with_headers=True — (data, None, headers_dict)
    с именами заголовков в нижнем регистре (например "x-total-results"), чтобы
    вызывающий мог прочитать служебные заголовки без забот о регистре.
    """
    backoff = 0.5
    for attempt in range(retries + 1):
        try:
            r = SESSION.get(url, params=params, timeout=TIMEOUT)
        except requests.exceptions.RequestException as e:
            # Сетевые сбои (обрыв, таймаут) — временные, стоит повторить.
            if attempt < retries:
                time.sleep(backoff)
                backoff *= 2
                continue
            # ВАЖНО: при with_headers=True обязаны вернуть РОВНО три значения,
            # иначе вызывающий (search_uniprot_candidates) распаковывает
            # результат в три переменные и падает с ValueError.
            return (None, str(e), None) if with_headers else (None, str(e))

        if r.status_code in (429, 500, 502, 503, 504) and attempt < retries:
            time.sleep(backoff)
            backoff *= 2
            continue

        if 400 <= r.status_code < 500 and r.status_code != 429:
            # 4xx — постоянная ошибка, повторять бессмысленно. Чаще всего это
            # 400 от UniProt на кривые символы в запросе (скобки, кавычки):
            # по фразе «не принял запрос» вызывающий (fetch_uniprot,
            # search_uniprot_candidates) узнаёт этот случай и пробует
            # следующий вариант запроса. Три значения при with_headers=True —
            # см. комментарий выше.
            msg = (
                f"сервер не принял запрос (ошибка {r.status_code}) — обычно "
                "так бывает из-за скобок или кавычек в названии"
            )
            return (None, msg, None) if with_headers else (None, msg)

        try:
            r.raise_for_status()
        except requests.exceptions.RequestException as e:
            # Три значения при with_headers=True — см. комментарий выше.
            return (None, str(e), None) if with_headers else (None, str(e))

        try:
            data = r.json()
        except ValueError as e:
            msg = f"сервер вернул ответ, который не получилось разобрать как JSON: {e}"
            return (None, msg, None) if with_headers else (None, msg)
        return (data, None, {k.lower(): v for k, v in r.headers.items()}) if with_headers else (data, None)

    return (None, "не удалось выполнить запрос", None) if with_headers else (
        None,
        "не удалось выполнить запрос",
    )


# ---------------------------------------------------------------------------
# 1. UniProt
# ---------------------------------------------------------------------------

def fetch_uniprot(term):
    """Ищет белок в UniProt и возвращает запись-победителя.

    Поисковый эндпоинт UniProt возвращает сразу полную запись (функцию,
    последовательность, ссылки на PDB и т.д.) — отдельный запрос не нужен.

    Возвращает (dict, None) или (None, сообщение_об_ошибке).
    """
    log("[1/3] UniProt: ищу белок...")

    if UNIPROT_ACCESSION_RE.match(term):
        queries = [f"accession:{term}"]
    else:
        # Порядок запросов от точного к неточному:
        #   1) gene_exact — точное имя гена, лучший контроль результата;
        #   2) текстовый поиск по проверенным записям;
        #   3) текстовый поиск без ограничений.
        queries = [
            f"gene_exact:{term} AND reviewed:true",
            f"{term} AND reviewed:true",
            term,
        ]

    # Способ, которым найдена запись, определяет оценку уверенности:
    #   accession — точный UniProt-код;
    #   gene      — совпадение с ОСНОВНЫМ именем гена;
    #   synonym   — совпадение только с синонимом;
    #   text      — обычный текстовый поиск (ген не проверялся).
    match_mode = None
    human_priority_used = False
    data = None
    last_err = None
    for q in queries:
        data, err = safe_get(
            "https://rest.uniprot.org/uniprotkb/search",
            # 10 записей вместо 1: показываем неоднозначность запроса
            # (инсулин есть у мыши, крысы, человека) и есть из чего
            # отбирать по точному имени гена.
            params={"query": q, "format": "json", "size": 10},
        )
        if err:
            last_err = err
            if "не принял запрос" in err:
                # 400-класс: кривой сам запрос (например «p53)» со скобкой) —
                # это не сбой сети, пробуем следующий, более простой вариант.
                data = None
                continue
            return None, f"UniProt: {err}"
        results = data.get("results") or []
        if not results:
            continue
        if q.startswith("gene_exact:"):
            # Приоритет кандидатов (важно для синонимов: у человека HER2 —
            # синоним гена ERBB2, а у других организмов HER2 может быть
            # главным именем гена):
            #   1) человек, совпадение с основным именем гена;
            #   2) человек, совпадение с синонимом;
            #   3) любой организм, основное имя гена;
            #   4) любой организм по синониму.
            primary = [r for r in results if gene_primary_matches(r, term)]
            by_synonym = [r for r in results if gene_synonym_matches(r, term)]
            human_primary = [r for r in primary if is_human(r)]
            human_synonyms = [r for r in by_synonym if is_human(r)]
            candidates = (human_primary or human_synonyms or primary
                          or by_synonym)
            if not candidates:
                continue
            chosen = candidates[0]
            # Режим определяем по факту выбранной записи: выбрали запись,
            # где запрос — синоним, значит уверенность ниже.
            match_mode = "gene" if gene_primary_matches(chosen, term) else "synonym"
            # Приоритет человека реально повлиял на выбор, только если
            # человек выбран из нескольких кандидатов, где были и не-люди.
            pool = primary + by_synonym
            pool_human = [r for r in pool if is_human(r)]
            human_priority_used = (
                is_human(chosen) and len(pool) > 1
                and 0 < len(pool_human) < len(pool)
            )
            chosen_ids = {id(r) for r in candidates}
            rest = [r for r in results if id(r) not in chosen_ids]
            data["results"] = candidates + rest
        else:
            match_mode = "accession" if q.startswith("accession:") else "text"
        break

    if not data or not data.get("results"):
        if data is None and last_err and "не принял запрос" in last_err:
            # Все варианты запроса отвергнуты сервером (400) — почти всегда
            # это скобки или кавычки в названии. Говорим по-человечески.
            return None, (
                f'UniProt не принял запрос «{term}» — обычно так бывает из-за '
                "скобок, кавычек или других знаков в названии. Проверь написание."
            )
        return None, f'В UniProt не нашлось ничего по запросу "{term}"'

    results = data["results"]
    entry = results[0]
    other_matches = [
        f"{r.get('primaryAccession')} ({(r.get('organism') or {}).get('scientificName', '?')})"
        for r in results[1:]
    ]

    # UniProt может прислать значение null вместо отсутствия ключа — .get(key, {})
    # в таком случае вернёт None, и следующий .get() упадёт. Поэтому везде
    # "or {}" / "or []" после .get(), а не дефолт внутри него.
    protein_desc = entry.get("proteinDescription") or {}
    recommended_name = protein_desc.get("recommendedName") or {}
    # Элемент списка может сам быть null ([null] вместо []) — отфильтровываем.
    submission_names = [s for s in (protein_desc.get("submissionNames") or []) if s] or [{}]
    full_name = (
        (recommended_name.get("fullName") or {}).get("value")
        or (submission_names[0].get("fullName") or {}).get("value")
        or "(название не указано)"
    )

    gene_names = []
    for g in entry.get("genes") or []:
        if not g:
            continue
        value = (g.get("geneName") or {}).get("value")
        if value:
            gene_names.append(value)

    # Из текстов функции берём первый непустой (пустые словари пропускаем).
    function_texts = []
    for comment in entry.get("comments") or []:
        if not comment or comment.get("commentType") != "FUNCTION":
            continue
        for text in comment.get("texts") or []:
            value = (text or {}).get("value")
            if value:
                function_texts.append(value)
                break

    # PubMed-ссылки встречаются в текстах комментариев (функция, активность,
    # заболевания и т.д.) в виде "PubMed:11025664". Собираем их отовсюду,
    # без повторов, в порядке первого появления — это литература, которой
    # UniProt подкрепляет данные записи.
    pubmed_ids = []
    seen_pubmed = set()
    for comment in entry.get("comments") or []:
        if not comment:
            continue
        for text in comment.get("texts") or []:
            value = (text or {}).get("value") or ""
            for match in re.finditer(r"PubMed:(\d+)", value):
                pmid = match.group(1)
                if pmid not in seen_pubmed:
                    seen_pubmed.add(pmid)
                    pubmed_ids.append(pmid)

    organism = entry.get("organism") or {}
    sequence = entry.get("sequence") or {}
    # length может быть null — тогда надёжнее посчитать самим.
    sequence_value = sequence.get("value") or ""
    sequence_length = sequence.get("length")
    if sequence_length is None:
        sequence_length = len(sequence_value)

    pdb_ids = [
        xref["id"]
        for xref in (entry.get("uniProtKBCrossReferences") or [])
        if xref and xref.get("database") == "PDB" and xref.get("id")
    ]

    return {
        "источник": "UniProt",
        "ссылка_на_запись": f"https://www.uniprot.org/uniprotkb/{entry.get('primaryAccession')}",
        "accession": entry.get("primaryAccession"),
        "id_записи": entry.get("uniProtkbId"),
        "название_белка": full_name,
        "гены": gene_names,
        "организм": organism.get("scientificName"),
        "организм_бытовое_название": organism.get("commonName"),
        "taxid": organism.get("taxonId"),
        "функция": function_texts[0] if function_texts else "(описание функции не найдено)",
        "литература": [
            {
                "pubmed_id": pmid,
                "ссылка": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            }
            for pmid in pubmed_ids
        ],
        "длина_последовательности_ак": sequence_length,
        "начало_последовательности": sequence_value[:40] + ("..." if len(sequence_value) > 40 else ""),
        "статус_проверки": entry.get("entryType"),
        # Как запись была найдена — для блока "УВЕРЕННОСТЬ" в карточке.
        "способ_поиска": {
            "accession": "точный UniProt-код",
            "gene": "точное совпадение с основным именем гена",
            "synonym": "совпадение с синонимом гена",
            "text": "обычный текстовый поиск (имя гена не проверялось)",
        }.get(match_mode, "неизвестно"),
        # Сырой режим поиска ("accession"/"gene"/"synonym"/"text") — нужен
        # batch-таблице: она различает «точное совпадение» (остальные записи —
        # просто то же имя у других организмов) и настоящий риск сопоставления.
        "режим_поиска": match_mode,
        "человек_выбран_автоматически": human_priority_used,
        "число_pdb_структур": len(pdb_ids),
        "pdb_id_список": pdb_ids,
        "другие_похожие_записи": other_matches,
    }, None


def search_uniprot_candidates(term):
    """Быстрый поиск по UniProt: только список кандидатов, без NCBI и PDB.

    Нужен для интерфейса: если по запросу несколько похожих записей,
    программа сначала показывает список и даёт пользователю выбрать самому,
    вместо того чтобы молча брать первую запись. Научная честность: выбор
    записи — решение человека, а не программы.

    Порядок кандидатов тот же, что при автоматическом выборе внутри
    fetch_uniprot (человек по главному имени гена > человек по синониму >
    любой организм по главному имени > по синониму > остальные), поэтому
    список начинается с самых вероятных записей — но ни одна из них НЕ
    подставляется автоматически.

    Для каждого кандидата заполняется «способ_поиска» — как именно он
    нашёлся по ЗАПРОСУ человека. Интерфейс передаёт это в build_card после
    выбора пользователя: без этого карточка по выбранному коду выглядела бы
    для assess_confidence как «введён точный UniProt-код», и шкала
    уверенности всегда показывала бы «высокую». Факт ручного выбора из
    списка передаётся ОТДЕЛЬНЫМ флагом (user_confirmed) и способ не
    подменяет — иначе почти любой поиск получал бы «высокую».

    Возвращает (список_кандидатов, None, всего_найдено) или
    (None, сообщение_об_ошибке, None). всего_найдено — общее число подходящих
    записей по данным UniProt (заголовок x-total-results), может быть None.
    """
    if UNIPROT_ACCESSION_RE.match(term):
        queries = [f"accession:{term}"]
    else:
        queries = [
            f"gene_exact:{term} AND reviewed:true",
            f"{term} AND reviewed:true",
            term,
        ]

    total_found = None
    last_err = None
    for q in queries:
        data, err, headers = safe_get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={"query": q, "format": "json", "size": 10},
            with_headers=True,
        )
        if err:
            last_err = err
            if "не принял запрос" in err:
                # 400-класс — кривой сам запрос (скобки, кавычки): пробуем
                # следующий вариант, а не обрываем весь поиск.
                continue
            return None, f"UniProt: {err}", None
        results = data.get("results") or []
        if not results:
            continue

        if q.startswith("gene_exact:"):
            primary = [r for r in results if gene_primary_matches(r, term)]
            by_synonym = [r for r in results if gene_synonym_matches(r, term)]
            human_primary = [r for r in primary if is_human(r)]
            human_synonyms = [r for r in by_synonym if is_human(r)]
            ordered = (human_primary or human_synonyms or primary or by_synonym)
            ordered_ids = {id(r) for r in ordered}
            results = ordered + [r for r in results if id(r) not in ordered_ids]

        # Общее число подходящих записей: UniProt отдаёт его в заголовке.
        # size жёстко ограничен десятью, поэтому список — это только первые
        # десять, и честно говорить об этом в интерфейсе.
        try:
            total_found = int(headers.get("x-total-results", ""))
        except (TypeError, ValueError):
            total_found = None

        candidates = []
        for r in results:
            if not r:
                continue
            # Разбор названия — те же правила устойчивости к null, что и в
            # fetch_uniprot (UniProt может прислать null вместо отсутствия).
            protein_desc = r.get("proteinDescription") or {}
            recommended = (protein_desc.get("recommendedName") or {}).get("fullName") or {}
            submission_names = [s for s in (protein_desc.get("submissionNames") or []) if s] or [{}]
            full_name = (
                recommended.get("value")
                or (submission_names[0].get("fullName") or {}).get("value")
                or "(название не указано)"
            )
            genes = []
            for g in r.get("genes") or []:
                if not g:
                    continue
                value = (g.get("geneName") or {}).get("value")
                if value:
                    genes.append(value)
            organism = r.get("organism") or {}
            if q.startswith("accession:"):
                way = "точный UniProt-код"
            elif q.startswith("gene_exact:"):
                if gene_primary_matches(r, term):
                    way = "точное совпадение с основным именем гена"
                elif gene_synonym_matches(r, term):
                    way = "совпадение с синонимом гена"
                else:
                    way = "обычный текстовый поиск (имя гена не проверялось)"
            else:
                way = "обычный текстовый поиск (имя гена не проверялось)"
            candidates.append({
                "accession": r.get("primaryAccession"),
                "название_белка": full_name,
                "организм": organism.get("scientificName"),
                "организм_отображение": organism_display(organism.get("scientificName")),
                "гены": genes,
                "статус_проверки": r.get("entryType"),
                "способ_поиска": way,
            })
        return candidates, None, total_found

    # Дошли сюда — ни один вариант запроса не дал записей. Если при этом
    # сервер отверг сами запросы (400) — объясняем по-человечески.
    if last_err and "не принял запрос" in last_err:
        return None, (
            f'UniProt не принял запрос «{term}» — обычно так бывает из-за '
            "скобок, кавычек или других знаков в названии. Проверь написание."
        ), None
    return [], None, None


# ---------------------------------------------------------------------------
# 2. NCBI (Entrez eutils)
# ---------------------------------------------------------------------------

def fetch_ncbi_gene(gene_symbol, taxid):
    """Ищет ген в NCBI Gene по символу гена + taxid организма (из UniProt).

    ВАЖНО: NCBI Gene хранит гены только организмов с аннотированным геномом.
    Для многих белков (например GFP медузы) записи гена нет вовсе — это не
    ошибка запроса, а реальный пробел в данных, который связывающий слой
    должен честно показывать.
    """
    log("[2/3] NCBI Entrez (gene): ищу ген...")

    if not gene_symbol:
        return None, "Нет символа гена от UniProt, поиск в NCBI Gene пропущен"

    term = f"{gene_symbol}[sym]"
    if taxid:
        term += f" AND {taxid}[taxid]"

    data, err = safe_get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params=ncbi_params(db="gene", term=term, retmax=1),
    )
    if err:
        return None, f"NCBI: {err}"

    idlist = (data.get("esearchresult") or {}).get("idlist") or []
    if not idlist:
        return None, f'В NCBI Gene не нашлось записи для "{term}" (для этого организма ген может быть неаннотирован)'

    time.sleep(ncbi_delay())

    gene_id = idlist[0]
    data2, err2 = safe_get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        params=ncbi_params(db="gene", id=gene_id),
    )
    if err2:
        return None, f"NCBI: {err2}"

    g = (data2.get("result") or {}).get(gene_id) or {}
    if not g.get("name"):
        # esummary может ответить 200 OK, но без данных по gene_id.
        return None, f"NCBI Gene вернул пустой ответ для Gene ID {gene_id}"

    return {
        "источник": "NCBI Gene",
        "ссылка_на_запись": f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}",
        "gene_id": gene_id,
        "символ": g.get("name"),
        "описание": g.get("description"),
        "хромосома": g.get("chromosome"),
        "расположение_в_геноме": g.get("maplocation"),
        "другие_названия": g.get("otheraliases"),
        "краткое_описание": (g.get("summary") or "(нет описания)")[:500],
    }, None


def fetch_ncbi_protein_fallback(term, organism_name):
    """Запасной вариант, если записи гена в NCBI нет: поиск в NCBI Protein.

    Это не запись гена (нет хромосомы/расположения), но подтверждает наличие
    белка в базах NCBI и даёт GenBank ID.
    """
    log("      -> запись гена не найдена, пробую NCBI Protein как запасной вариант...")

    query = term if not organism_name else f"{term} {organism_name}"
    data, err = safe_get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params=ncbi_params(db="protein", term=query, retmax=1),
    )
    if err:
        log(f"      -> NCBI Protein: {err}")
        return None
    idlist = (data.get("esearchresult") or {}).get("idlist") or []
    if not idlist:
        log('      -> NCBI Protein: тоже ничего не нашлось')
        return None

    time.sleep(ncbi_delay())
    protein_id = idlist[0]
    data2, err2 = safe_get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        params=ncbi_params(db="protein", id=protein_id),
    )
    if err2:
        log(f"      -> NCBI Protein: {err2}")
        return None

    p = (data2.get("result") or {}).get(protein_id) or {}
    return {
        "источник": "NCBI Protein (запасной вариант, т.к. записи гена нет)",
        "ссылка_на_запись": f"https://www.ncbi.nlm.nih.gov/protein/{protein_id}",
        "genbank_id": p.get("accessionversion"),
        "заголовок": p.get("title"),
        "организм": p.get("organism"),
        "длина_ак": p.get("slen"),
        # Это текстовый поиск по названию без проверки "тот самый ли ген",
        # поэтому доверять этому блоку меньше, чем ncbi_gene.
        "надёжность": "низкая — первый результат текстового поиска, не проверен вручную",
    }


# ---------------------------------------------------------------------------
# 3. RCSB PDB — 3D-структуры
# ---------------------------------------------------------------------------

def fetch_pdb_structures(pdb_ids, max_structures=3):
    """Для первых PDB ID из UniProt запрашивает детали структур в RCSB."""
    log("[3/3] RCSB PDB: смотрю 3D-структуры...")

    if not pdb_ids:
        # Мы проверили не весь RCSB, а только то, что о нём знает UniProt
        # (cross-reference) — это не одно и то же.
        return None, "В UniProt не найдено ссылок на PDB-структуры для этой записи"

    structures = []
    fetch_errors = []
    for pdb_id in pdb_ids[:max_structures]:
        data, err = safe_get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
        if err:
            fetch_errors.append(f"{pdb_id}: {err}")
            continue
        exptl = data.get("exptl") or [{}]
        entry_info = data.get("rcsb_entry_info") or {}
        structures.append({
            "pdb_id": pdb_id,
            "ссылка": f"https://www.rcsb.org/structure/{pdb_id}",
            "название": (data.get("struct") or {}).get("title"),
            "метод": (exptl[0] or {}).get("method") if exptl else None,
            "разрешение_ангстрем": (entry_info.get("resolution_combined") or [None])[0],
        })

    if not structures:
        return None, "Не удалось получить детали структур из RCSB: " + "; ".join(fetch_errors)

    return {
        "источник": "RCSB PDB",
        "всего_структур_в_uniprot": len(pdb_ids),
        "показано_примеров": len(structures),
        "структуры": structures,
        "ошибки_запросов": fetch_errors,
    }, None


# ---------------------------------------------------------------------------
# Склейка карточки
# ---------------------------------------------------------------------------

def assess_confidence(uniprot, search_way=None, user_confirmed=False):
    """Оценивает надёжность сопоставления "запрос → найденная запись".

    Главная научная ценность подхода — проверка, что из трёх баз собраны
    сведения про ОДИН и ТОТ ЖЕ объект. Карточка честно сообщает, насколько
    мы уверены, что нашли именно то, что искали.

    search_way — способ поиска, известный интерфейсу ДО сборки карточки
    (когда человек сам выбрал запись из списка кандидатов). Без него оценка
    считала бы способ по факту: «карточка собрана по точному коду UniProt» —
    и всегда возвращала «высокая», даже если запись нашлась слабым текстовым
    поиском.

    user_confirmed — человек сам выбрал запись из списка похожих. Это
    ОТДЕЛЬНОЕ измерение, а не уровень шкалы: выбор человека не подменяет
    способ сопоставления. Кандидатов больше одного почти всегда (у insulin
    их 5435), и если бы выбор засчитывался как «высокая», почти любой поиск
    получал бы зелёную плашку — и шкала переставала бы что-либо означать.
    Поэтому факт ручного выбора честно дописывается к причине, а уровень
    определяет только способ, которым запись нашлась.
    """
    way = search_way or uniprot.get("способ_поиска", "")
    human_used = uniprot.get("человек_выбран_автоматически")

    if way.startswith("точный UniProt-код"):
        level, reason = "высокая", "введён точный UniProt-код — запись определена однозначно."
    elif way.startswith("точное совпадение"):
        if human_used:
            level, reason = "высокая", (
                "запрос совпал с основным именем гена; организмов с таким именем"
                " несколько — выбран организм человека (остальные перечислены в"
                " похожих записях)"
            )
        else:
            level, reason = "высокая", "запрос совпал с основным именем гена."
    elif way.startswith("совпадение с синонимом"):
        if human_used:
            level, reason = "средняя", (
                "запрос — синоним гена; из нескольких организмов выбран человек"
                " (автоматический приоритет)"
            )
        else:
            level, reason = "средняя", "запрос — синоним гена, возможны другие организмы."
    else:
        # Текстовый поиск (или неизвестный способ) — самый слабый случай:
        # релевантность решил сам UniProt. Здесь ручной выбор особенно важен
        # упомянуть: без него организм выбирал бы поиск, а не человек.
        if user_confirmed:
            return "низкая", (
                "запись найдена обычным текстовым поиском; совпадение имени гена"
                " не проверялось. Человек сам выбрал её из списка похожих —"
                " организм подтверждён вручную."
            )
        return "низкая", (
            "запись найдена обычным текстовым поиском; совпадение имени гена"
            " не проверялось — какой организм показать, решил сам поиск."
        )

    if user_confirmed:
        reason = reason.rstrip(".") + (
            ". Выбор записи подтверждён вручную — человек выбрал её из списка"
            " похожих."
        )
    return level, reason


def build_card(term, search_way=None, user_confirmed=False, детали_pdb=True):
    """Собирает карточку: UniProt → NCBI → PDB, с честными ошибками.

    search_way, user_confirmed — см. assess_confidence: сайт передаёт
    способ, которым выбранная человеком запись нашлась в списке кандидатов,
    и флаг «человек выбрал запись сам». Флаг уровень не подменяет — только
    дополняет причину фактом ручного подтверждения.

    детали_pdb=False — быстрый режим массовой обработки: детали структур
    (название, метод, разрешение) в RCSB не запрашиваются — минус до трёх
    HTTP-запросов на каждый белок, примерно вдвое быстрее. Число структур
    всё равно приходит из UniProt, а детали в таблицу batch не попадают,
    поэтому для неё ничего не теряется.
    """
    card = {
        "запрос_пользователя": term,
        "время_сборки": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "источники": {},
        "проблемы_источников": [],
    }

    uniprot, err = fetch_uniprot(term)
    if uniprot:
        card["источники"]["uniprot"] = uniprot
        level, reason = assess_confidence(uniprot, search_way, user_confirmed)
        card["уверенность"] = {"уровень": level, "причина": reason}
    else:
        card["проблемы_источников"].append(err)
        log(f"      -> {err}")
        # Без UniProt нет ни гена, ни организма для дальнейшего поиска.
        return card

    # Пауза перед ПЕРВЫМ запросом к NCBI в этой карточке: предыдущая карточка
    # (например, в массовом режиме) могла закончить запросы к eutils только
    # что. Одна строчка делает темп вежливым и без API-ключа (~2 запроса/сек
    # при лимите 3) при любом порядке источников.
    time.sleep(ncbi_delay())
    gene_symbol = uniprot["гены"][0] if uniprot["гены"] else term
    ncbi_gene, ncbi_err = fetch_ncbi_gene(gene_symbol, uniprot.get("taxid"))
    if ncbi_gene:
        card["источники"]["ncbi_gene"] = ncbi_gene
    else:
        card["проблемы_источников"].append(ncbi_err)
        log(f"      -> {ncbi_err}")
        # Перед запасным запросом к тому же серверу NCBI нужна пауза.
        time.sleep(ncbi_delay())
        fallback = fetch_ncbi_protein_fallback(gene_symbol, uniprot.get("организм"))
        if fallback:
            card["источники"]["ncbi_protein_fallback"] = fallback

    if детали_pdb:
        pdb, pdb_err = fetch_pdb_structures(uniprot.get("pdb_id_список", []))
        if pdb:
            card["источники"]["pdb"] = pdb
        else:
            card["проблемы_источников"].append(pdb_err)
            log(f"      -> {pdb_err}")

    return card


# ---------------------------------------------------------------------------
# Массовая обработка (batch)
# ---------------------------------------------------------------------------

CSV_HEADER_WORDS = {
    "gene", "genes", "gen", "ген", "гена", "name", "names", "название", "названия",
    "имя", "query", "запрос", "protein", "proteins", "белок", "белки",
    "id", "ids", "identifier", "identifiers", "uniprot", "accession",
    "entry", "symbol", "from", "to", "код",
}


def looks_like_csv_header(line):
    """Первая строка csv-файла похожа на заголовок, если ВСЕ её слова —
    типичные заголовочные («Gene name», «Protein ID», «From»).

    Проверять только первое слово нельзя: настоящее название белка может на
    него начинаться («Protein kinase...»). А целиком из заголовочных слов
    настоящее название состоять не может. Цифры оставляем при слове:
    иначе настоящий ген ID2 разобрался бы на "id" + "2" и был бы принят
    за заголовок.
    """
    words = re.findall(r"[a-zа-яё0-9]+", line.lower())
    return bool(words) and all(w in CSV_HEADER_WORDS for w in words)


def parse_batch_input(raw, csv_mode=False):
    """Разбирает вставленный список или содержимое файла в список запросов.

    Разделители — перевод строки, запятая, точка с запятой. Пробел и таб
    разделителями НЕ считаются: названия белков бывают многословными
    («Cellular tumor antigen p53»), и резать их по пробелу нельзя.
    В csv_mode (загруженный .csv) из каждой строки берётся только первая
    колонка — иначе лишние колонки («human», «Organism») и строка заголовка
    превращаются в мусорные запросы; строка-заголовок («Gene name»,
    «Protein ID», «From») распознаётся по словам и выбрасывается. Поля в
    кавычках (Excel так оформляет значения с запятой) читаются модулем csv:
    название вроде "Kinase, alpha" остаётся целым.
    Точка в конце названия не отрезается: существуют версионные коды записей
    вроде P01308.1, и точку у них отрезать нельзя.
    Повторы убираются (регистронезависимо), порядок сохраняется.

    Возвращает (список_запросов, число_убранных_повторов).
    """
    if not raw:
        return [], 0
    terms = []
    seen = set()
    duplicates = 0
    for line in raw.strip().splitlines():
        if csv_mode:
            # Excel берёт поля с запятой в кавычки: «"Kinase, alpha",human».
            # re.split разрезал бы название по запятой внутри кавычек;
            # csv.reader понимает кавычки и отдаёт первую колонку целой.
            # Первая колонка csv — уже готовый термин, второй раз по
            # запятым её резать нельзя.
            try:
                first_row = next(csv.reader([line], skipinitialspace=True))
            except (csv.Error, StopIteration):
                first_row = []
            # Берём ровно первую колонку: остальное («human», «Organism») —
            # мусор для запроса. Пустая строка → обрабатываем как есть.
            pieces = [first_row[0]] if first_row else [line]
        else:
            pieces = re.split(r"[,;]+", line)
        for piece in pieces:
            term = piece.strip().strip("\"'").strip()
            if not term:
                continue
            if csv_mode and not terms and looks_like_csv_header(term):
                continue  # строка-заголовок csv («Gene name», «Protein ID», «From»)
            key = term.lower()
            if key in seen:
                duplicates += 1
                continue
            seen.add(key)
            terms.append(term)
    return terms, duplicates


def batch_row_from_card(card):
    """Превращает карточку в плоскую строку для таблицы и CSV.

    Никогда не бросает исключений: если сборка карточки не удалась совсем,
    строка честно помечается "ошибка" с причиной.
    """
    term = card.get("запрос_пользователя", "")
    sources = card.get("источники") or {}
    up = sources.get("uniprot")
    conf = card.get("уверенность") or {}
    ncbi_gene = sources.get("ncbi_gene") or {}
    problems = card.get("проблемы_источников") or []

    row = {
        "запрос": term,
        "результат": "не найдено",
        "uniprot_accession": "",
        "название_белка": "",
        "организм": "",
        "гены": "",
        "уверенность": "",
        "статус_проверки": "",
        "хромосома_ncbi": "",
        "структур_в_pdb": "",
        "ссылка_uniprot": "",
        "примечание": "",
    }

    if up:
        row["результат"] = "найдено"
        row["uniprot_accession"] = up.get("accession", "")
        row["название_белка"] = up.get("название_белка", "")
        row["организм"] = organism_display(up.get("организм"))
        row["гены"] = ", ".join(up.get("гены") or [])
        row["уверенность"] = (conf.get("уровень") or "").upper()
        row["статус_проверки"] = entry_status_display(up.get("статус_проверки"))
        row["хромосома_ncbi"] = ncbi_gene.get("хромосома") or ""
        row["структур_в_pdb"] = up.get("число_pdb_структур", "")
        row["ссылка_uniprot"] = up.get("ссылка_на_запись", "")

        others = up.get("другие_похожие_записи") or []
        if others:
            shown = "; ".join(others[:3])
            more = f" ... (+{len(others) - 3})" if len(others) > 3 else ""
            organism = row["организм"] or "?"
            mode = up.get("режим_поиска")
            if mode in ("accession", "gene"):
                # Совпали точный код или основное имя гена: остальные записи —
                # просто то же имя у других организмов (у TP53 это p53 крысы,
                # мыши, рыбки данио…). Выбор человека среди организмов — не
                # неоднозначность, а нормальная работа приоритета, поэтому
                # пометка спокойная. Слово «неоднозначно» тут нарочно
                # отсутствует: summarize_batch считает спорные строки именно
                # по этому слову в примечании.
                what = (
                    "совпал точный код UniProt" if mode == "accession"
                    else "совпало основное имя гена"
                )
                row["примечание"] = (
                    f"взята первая — {organism}; {what}; такое имя есть ещё "
                    f"у: {shown}{more}"
                )
            else:
                # Совпадение по синониму или текстовым поиском — настоящий
                # риск сопоставления, честно помечаем как неоднозначно.
                how = (
                    "совпадение по синониму гена" if mode == "synonym"
                    else "текстовый поиск, имя гена не проверялось"
                )
                row["примечание"] = (
                    f"неоднозначно: взята первая — {organism}; {how}; есть "
                    f"ещё похожие записи: {shown}{more}"
                )
        elif problems:
            # Запись UniProt есть, но NCBI/PDB не ответили — это не ошибка
            # сопоставления, а пробелы в данных источников; показываем честно.
            row["примечание"] = "частично: " + "; ".join(problems[:2])
    else:
        note = problems[0] if problems else ""
        if note and "не нашлось" not in note:
            row["результат"] = "ошибка"
        row["примечание"] = note

    return row


def summarize_batch(rows):
    """Сводка по строкам batch: сколько найдено, по уровням уверенности и т.д."""
    summary = {
        "всего": len(rows),
        "найдено": 0,
        "высокая": 0,
        "средняя": 0,
        "низкая": 0,
        "не найдено": 0,
        "ошибок": 0,
        "неоднозначных": 0,
    }
    for row in rows:
        result = row.get("результат")
        if result == "найдено":
            summary["найдено"] += 1
            level = (row.get("уверенность") or "").lower()
            if level in summary:
                summary[level] += 1
            if "неоднозначно" in (row.get("примечание") or ""):
                summary["неоднозначных"] += 1
        elif result == "не найдено":
            summary["не найдено"] += 1
        elif result == "ошибка":
            summary["ошибок"] += 1
    return summary


# ---------------------------------------------------------------------------
# Вывод
# ---------------------------------------------------------------------------

def render_text(card):
    lines = []
    lines.append("=" * 70)
    lines.append(f'КАРТОЧКА ОБЪЕКТА: "{card["запрос_пользователя"]}"')
    lines.append(f'Собрано: {card["время_сборки"]}')
    conf = card.get("уверенность")
    if conf:
        lines.append(f'УВЕРЕННОСТЬ: {conf["уровень"].upper()} — {conf["причина"]}')
    lines.append("=" * 70)

    up = card["источники"].get("uniprot")
    if up:
        lines.append("")
        lines.append("--- UniProt (белок) ---")
        lines.append(f"Название: {up['название_белка']}")
        lines.append(f"Accession: {up['accession']} ({up['id_записи']})")
        lines.append(f"Статус проверки: {entry_status_display(up.get('статус_проверки'))}")
        lines.append(f"Гены: {', '.join(g for g in up['гены'] if g) or '(не указаны)'}")
        lines.append(f"Организм: {organism_display(up.get('организм'))}")
        lines.append(f"Функция: {up['функция']}")
        lines.append(f"Длина последовательности: {up['длина_последовательности_ак']} аминокислот")
        lines.append(f"Начало последовательности: {up['начало_последовательности']}")
        lines.append(f"Известных 3D-структур в PDB (по данным UniProt): {up['число_pdb_структур']}")
        lines.append(f"Ссылка: {up['ссылка_на_запись']}")
        others = up.get("другие_похожие_записи") or []
        if others:
            mode = up.get("режим_поиска")
            org = organism_display(up.get("организм"))
            if mode in ("accession", "gene"):
                lines.append(
                    f"i То же имя есть ещё у других организмов — выбран {org} "
                    "по приоритету (это не неоднозначность): "
                    + ", ".join(others)
                )
            else:
                lines.append(
                    f"! Запрос неоднозначный — взята первая ({org}), есть ещё "
                    "похожие записи: " + ", ".join(others)
                )
        lit = up.get("литература") or []
        if lit:
            shown = lit[:10]
            lines.append("")
            lines.append(f"Литература (PubMed): показано {len(shown)} из {len(lit)}")
            for pub in shown:
                lines.append(f"  * PMID {pub['pubmed_id']} — {pub['ссылка']}")
            if len(lit) > len(shown):
                lines.append("  ... полный список — в записи UniProt по ссылке выше")

    ng = card["источники"].get("ncbi_gene")
    if ng:
        lines.append("")
        lines.append("--- NCBI Gene (ген) ---")
        lines.append(f"Символ гена: {ng['символ']} (Gene ID {ng['gene_id']})")
        lines.append(f"Описание: {ng['описание'] or '(нет описания)'}")
        lines.append(f"Хромосома: {ng['хромосома'] or '—'} ({ng['расположение_в_геноме'] or '—'})")
        lines.append(f"Другие названия: {ng['другие_названия'] or '—'}")
        lines.append(f"Краткое описание: {ng['краткое_описание']}")
        lines.append(f"Ссылка: {ng['ссылка_на_запись']}")

    nf = card["источники"].get("ncbi_protein_fallback")
    if nf:
        lines.append("")
        lines.append("--- NCBI Protein (запасной вариант) ---")
        lines.append(f"GenBank ID: {nf['genbank_id']}")
        lines.append(f"Заголовок: {nf['заголовок']}")
        lines.append(f"Организм: {nf['организм']}")
        lines.append(f"Надёжность: {nf['надёжность']}")
        lines.append(f"Ссылка: {nf['ссылка_на_запись']}")

    pdb = card["источники"].get("pdb")
    if pdb:
        lines.append("")
        lines.append("--- RCSB PDB (3D-структуры) ---")
        lines.append(f"Всего структур по данным UniProt: {pdb['всего_структур_в_uniprot']}, показаны примеры: {pdb['показано_примеров']}")
        for s in pdb["структуры"]:
            # У NMR-структур разрешения не бывает (оно есть только у
            # рентгена и крио-ЭМ) — честно пишем "не применимо".
            resolution = s['разрешение_ангстрем'] if s['разрешение_ангстрем'] is not None else 'нет (не применимо)'
            lines.append(f"  * {s['pdb_id']}: {s['название'] or '(без названия)'} | метод: {s['метод'] or '—'} | разрешение: {resolution} | {s['ссылка']}")
        if pdb.get("ошибки_запросов"):
            lines.append(f"  (не удалось загрузить: {'; '.join(pdb['ошибки_запросов'])})")

    if card["проблемы_источников"]:
        lines.append("")
        lines.append("--- Что не удалось получить ---")
        for p in card["проблемы_источников"]:
            lines.append(f"  ! {p}")

    lines.append("")
    lines.append("=" * 70)
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        prog="prototype.py",
        description="Собирает карточку белка/гена из UniProt, NCBI Gene/Protein и RCSB PDB.",
        epilog=(
            "Примеры:\n"
            "  python3 prototype.py insulin\n"
            "  python3 prototype.py P01308\n"
            "  python3 prototype.py GFP"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "query",
        help="имя гена/белка (например insulin, GFP) или UniProt-код (например P01308)",
    )
    args = parser.parse_args()

    term = args.query.strip()
    if not term:
        print("Ошибка: пустой запрос. Укажи имя гена/белка или UniProt-код.")
        sys.exit(1)

    log(f'Собираю карточку по запросу: "{term}"\n')

    card = build_card(term)

    text = render_text(card)
    print("\n" + text)

    # Имя файла — по UniProt accession (короткий и однозначный). Если
    # accession нет (UniProt ничего не нашёл), используем очищенный запрос
    # с ограничением длины и защитой от пустого имени.
    accession = card.get("источники", {}).get("uniprot", {}).get("accession")
    safe_name = accession or re.sub(r"[^a-zA-Zа-яА-Я0-9_-]+", "_", term)[:MAX_SAFE_NAME_LEN]
    safe_name = safe_name.strip("_") or "запрос"
    txt_path = f"card_{safe_name}.txt"
    json_path = f"card_{safe_name}.json"

    with open(txt_path, "w", encoding="utf-8") as f:
        f.write(text + "\n")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(card, f, ensure_ascii=False, indent=2)

    print(f"\nКарточка сохранена в файлы:\n  {txt_path}\n  {json_path}")


if __name__ == "__main__":
    main()

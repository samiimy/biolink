# -*- coding: utf-8 -*-
"""
Веб-версия BioLink: собирает карточку белка/гена из UniProt, NCBI и PDB
и показывает её в браузере. Вся логика поиска живёт в
prototype_semifinale.py — этот файл только рисует интерфейс.

Две вкладки:
  1. "Поиск карточки" — один запрос -> карточка. Если по запросу
     находится несколько похожих записей, программа НЕ берёт первую молча,
     а сначала показывает список кандидатов и даёт пользователю выбрать.
  2. "Массовая обработка" — список названий (вставить или загрузить файл)
     -> программа по очереди собирает карточки и складывает всё в таблицу
     с уровнем уверенности; результат скачивается в CSV.

Запуск:  streamlit run app.py
"""

import json
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests
import streamlit as st

# Чтобы app.py можно было запускать из любой папки.
sys.path.insert(0, str(Path(__file__).parent))
import prototype_semifinale as core  # noqa: E402

# 3D-модель в браузере: py3Dmol собирает HTML-виджет на движке 3Dmol.js,
# показываем его через st.components.v1.html (это ровно то, что делает
# обёртка stmol — components.html(view._make_html()) — но без её тяжелих
# jupyter-зависимостей, которые конфликтуют с новыми версиями IPython).
# Пакет ставится отдельно (pip install py3Dmol), поэтому импорт мягкий:
# если его нет, сайт работает как раньше, а в карточке вместо модели —
# подсказка, как его установить.
try:
    import py3Dmol
    RENDER3D_OK = True
except ImportError:
    RENDER3D_OK = False


def load_ncbi_secrets():
    """API-ключ NCBI из .streamlit/secrets.toml (файл не лежит в коде).

    Без ключа NCBI разрешает ~3 запроса/сек, с бесплатным ключом — до 10.
    """
    try:
        key = st.secrets["ncbi"]["api_key"]
        if key:
            core.NCBI_API_KEY = str(key).strip()
    except Exception:
        pass


load_ncbi_secrets()


def linkify_pubmed(text):
    """Делает ссылки PubMed:XXXX в текстах UniProt кликабельными."""
    return re.sub(
        r"PubMed:(\d+)",
        r"[PubMed:\1](https://pubmed.ncbi.nlm.nih.gov/\1/)",
        text,
    )


def format_duration(seconds):
    """Человекочитаемое время: '4 мин' / '1 ч 20 мин'."""
    seconds = int(seconds)
    if seconds < 60:
        return f"~{seconds} сек"
    minutes = seconds // 60
    if minutes < 60:
        return f"~{minutes} мин"
    hours, minutes = divmod(minutes, 60)
    return f"~{hours} ч {minutes} мин"


# Средняя длительность сборки одной карточки в массовом режиме
# (ядро в этом режиме не запрашивает детали структур в PDB: 2–3 секунды).
EST_SECONDS_PER_CARD = 3


@st.cache_data(ttl=24 * 3600, show_spinner=False, max_entries=32)
def cached_pdb_text(pdb_id: str) -> str:
    """Файл структуры в формате PDB с сервера RCSB; кэш на сутки.

    Скачиваем сами (а не через query= в py3Dmol), чтобы повторные показы
    модели не дёргали сеть: один и тот же файл подходит и для вращения,
    и для повторного включения тумблера.

    max_entries=32: файл структуры весит до нескольких МБ, и без лимита
    кэш разрастался бы после каждого просмотра десятка структур.
    """
    resp = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=30)
    resp.raise_for_status()
    return resp.text


def render_3d_block(pdb, up):
    """Интерактивная 3D-модель белка: вращение мышкой (3Dmol.js + py3Dmol).

    Спрятана за тумблером: пока его не включили, ничего не скачивается и
    не рендерится — открытие карточки не становится медленнее. Рендер идёт
    только на вкладке «Поиск»: в массовом режиме скорость важнее.

    accession в ключах виджетов ОБЯЗАТЕЛЕН: если ключ постоянный, Streamlit
    при переходе к другому белку подсовывает выбранную ранее структуру
    (её нет в новом списке, но старое значение возвращается молча) —
    и на экране крутится чужой белок, подписанный как текущий.
    """
    with st.expander("🧬 3D-модель — покрутить мышкой"):
        if not RENDER3D_OK:
            st.info(
                "Модель не показывается: не установлен пакет py3Dmol. "
                "Команда: `pip install py3Dmol`"
            )
            return

        accession = up["accession"]
        # Выбирать можно из ВСЕХ структур карточки UniProt, а не только из
        # первых с запрошенными деталями: подписи (название) есть у детальных,
        # остальные показываются просто кодом — файл структуры всё равно
        # качается по коду с files.rcsb.org, детали для рендера не нужны.
        all_ids = up.get("pdb_id_список") or []
        detailed = (pdb or {}).get("структуры") or []
        labels = {
            s["pdb_id"]: f"{s['pdb_id']} — {s['название'] or '(без названия)'}"
            for s in detailed
        }
        col_a, col_b = st.columns([3, 2])
        pdb_id = col_a.selectbox(
            f"Структура (всего {len(all_ids)})",
            list(all_ids),
            format_func=lambda p: labels.get(p, p),
            key=f"pdb3d_структура_{accession}",
        )
        style_name = col_b.selectbox(
            "Стиль",
            ["лента (cartoon)", "палочки (stick)", "шары (sphere)", "линии (line)"],
            key=f"pdb3d_стиль_{accession}",
        )

        if not st.toggle("Показать 3D-модель", key=f"pdb3d_показать_{accession}"):
            st.caption(
                "Включи тумблер — структура загрузится с сервера RCSB и появится здесь. "
                "Вращение: зажать левую кнопку мыши; приблизить: колесо."
            )
            return

        try:
            with st.spinner("Загружаю структуру с RCSB…"):
                pdb_text = cached_pdb_text(pdb_id)
        except Exception:
            st.warning(
                f"Не удалось скачать структуру {pdb_id}. Открой её прямо на сайте: "
                f"[rcsb.org/structure/{pdb_id}](https://www.rcsb.org/structure/{pdb_id})."
            )
            return

        # Крио-ЭМ структуры бывают по десяткам мегабайт — такой файл
        # тормозит даже хороший браузер, честнее отправить на сайт RCSB.
        if len(pdb_text) > 8_000_000:
            st.warning(
                "Эта структура очень большая и может тормозить в браузере. "
                f"Если модель не крутится, открой её на сайте: "
                f"[rcsb.org/structure/{pdb_id}](https://www.rcsb.org/structure/{pdb_id})."
            )

        # Из подписи «лента (cartoon)» берём ключ стиля для 3Dmol.js.
        style_key = style_name.split("(")[1].rstrip(")")
        # Ширина '100%' (строка, py3Dmol вставляет её в CSS как есть): холст
        # растягивается на всю ширину колонки и НЕ вылезает на телефонах.
        # Фиксированные пиксели (640) ломали мобильную вёрстку — модель
        # уезжала за правый край узкого экрана.
        view = py3Dmol.view(width="100%", height=480)
        view.addModel(pdb_text, "pdb")
        view.setStyle({style_key: {"color": "spectrum"}})
        view.zoomTo()
        # Тот же показ, что делает stmol.showmol: HTML-виджет 3Dmol.js в iframe.
        # iframe растянут на всю ширину карточки; холст 3Dmol внутри — тоже.
        st.components.v1.html(view._make_html(), height=490)
        st.caption(
            "Движок: [3Dmol.js](https://3dmol.csb.pitt.edu). "
            "Данные структуры: RCSB PDB."
        )

st.set_page_config(
    page_title="BioLink: карточка белка",
    layout="wide",
)

st.title("BioLink")
st.caption(
    "Вводишь имя белка, гена или код UniProt — программа спрашивает три научные базы "
    "(UniProt, NCBI, RCSB PDB) и собирает из них одну карточку. На вкладке "
    "«Массовая обработка» можно прогнать сразу список названий."
)


# ---------------------------------------------------------------------------
# Вкладка 1. Поиск одной карточки
# ---------------------------------------------------------------------------

def render_search_tab():
    EXAMPLES = ["insulin", "TP53", "p53", "HER2", "GFP", "P01308"]

    # Кнопки примеров — ДО поля ввода: клик пишет значение в состояние
    # заранее, пока поле ещё не создано (иначе Streamlit не даст менять
    # уже нарисованное поле), и поиск запускается сразу.
    st.write("**Быстрые примеры:**")
    cols = st.columns(len(EXAMPLES))
    clicked_example = None
    for col, example in zip(cols, EXAMPLES):
        if col.button(example, use_container_width=True):
            st.session_state["запрос"] = example
            clicked_example = example

    # Поле и кнопка внутри st.form: тогда Enter в поле запускает поиск так
    # же, как кнопка. Без формы Enter просто перерисовывал страницу,
    # и поиск не запускался.
    with st.form("поиск_карточки", border=False):
        query = st.text_input(
            "Что ищем?",
            placeholder="например: insulin, TP53, p53, GFP или код P01308",
            key="запрос",
        )
        run = st.form_submit_button("Собрать карточку", type="primary")

    term = clicked_example or (query or "").strip()

    # Кнопка в Streamlit "нажата" только в том прогоне, в котором по ней был
    # клик. Если пользователь выбрал запись из списка на Шаге выбора, страница
    # перерисуется уже без нажатой кнопки — поэтому активный запрос хранится
    # в session_state: выбор из списка перерисовывает карточку сам, без
    # повторного клика.
    if term and (clicked_example or run):
        st.session_state["активный_запрос"] = term
        if st.session_state.get("термин_выбора") != term:
            st.session_state["термин_выбора"] = term
            st.session_state.pop("выбор_записи", None)

    active_term = st.session_state.get("активный_запрос")
    if not active_term:
        if run and not term:
            st.error("Пустой запрос — введи имя гена/белка или код UniProt.")
        else:
            # Пустое состояние: без подсказки сайт при первом заходе молчит.
            st.info("Введи название белка, гена или код UniProt — или нажми на один из примеров выше 👆")
        return
    term = active_term

    @st.cache_data(ttl=6 * 3600, show_spinner=False)
    def _candidates_cached(query_text: str):
        """Только поиск кандидатов в UniProt (быстро, без NCBI и PDB).

        Ошибка превращается в исключение: st.cache_data НЕ сохраняет
        поднятое исключение, поэтому сетевой сбой не «консервируется» на
        все 6 часов ttl — следующий запрос повторит попытку по-настоящему.
        """
        candidates, err, total_found = core.search_uniprot_candidates(query_text)
        if err:
            raise RuntimeError(err)
        return candidates, total_found

    def cached_candidates(query_text: str):
        """Обёртка: возвращает обычный кортеж (кандидаты, ошибка, всего)."""
        try:
            candidates, total_found = _candidates_cached(query_text)
            return candidates, None, total_found
        except Exception as e:
            return None, str(e), None

    @st.cache_data(ttl=6 * 3600, show_spinner=False)
    def _card_cached(uniprot_accession: str, search_way: str, user_confirmed: bool):
        """Полная карточка; кэш на 6 часов — повторные запросы не дёргают базы.

        search_way — как ЭТА запись нашлась по запросу человека на шаге
        поиска (кандидат помнит свой способ), user_confirmed — выбрал ли
        человек её сам из списка. Оба входят в ключ кэша, чтобы карточки
        разных путей не смешивались. Без search_way оценка уверенности
        считала бы по факту «собрано по точному коду» и всегда писала
        «высокую» — см. core.assess_confidence.

        Карточка вообще без источников (UniProt не ответил) — это ошибка,
        а не результат: поднимаем исключение, чтобы она не осталась в кэше
        на все 6 часов. Частичные сбои (UniProt есть, NCBI/PDB нет) —
        обычный результат, кэшируются как есть.
        """
        card = core.build_card(
            uniprot_accession, search_way=search_way, user_confirmed=user_confirmed
        )
        if not card.get("источники"):
            raise RuntimeError(
                (card.get("проблемы_источников") or ["источники не ответили"])[0]
            )
        return card

    def cached_card(uniprot_accession: str, search_way: str = "", user_confirmed: bool = False):
        """Обёртка: ошибки возвращает карточкой без источников (не кэшируя)."""
        try:
            return _card_cached(uniprot_accession, search_way, user_confirmed)
        except Exception as e:
            return {
                "запрос_пользователя": uniprot_accession,
                "источники": {},
                "проблемы_источников": [str(e)],
            }

    with st.spinner("Ищу в UniProt…"):
        candidates, err, total_found = cached_candidates(term)

    if err:
        st.error(f"Ошибка обращения к UniProt: {err}")
        return
    if not candidates:
        st.error(
            "Ничего не нашлось в UniProt по этому запросу. "
            "Попробуй другое имя гена или UniProt-код."
        )
        return

    # ----------------------------------------------------------------------
    # Шаг выбора. Если записей несколько — программа не решает молча за
    # пользователя: сначала список кандидатов, потом карточка по выбранной.
    # Если запись одна (или введён точный код UniProt) — лишний шаг не нужен.
    # ----------------------------------------------------------------------
    # Как каждый кандидат был найден по запросу — нужно для честной шкалы
    # уверенности: карточка собирается по коду, и без этой подсказки сайт
    # считал бы любой выбор «точным кодом» (всегда «высокая»).
    way_by_accession = {c["accession"]: c.get("способ_поиска", "") for c in candidates}

    accession = candidates[0]["accession"]
    if len(candidates) > 1:
        total_note = ""
        if total_found and total_found > len(candidates):
            total_note = (
                f" Всего подходящих записей: {total_found} — показаны первые "
                f"{len(candidates)}; чтобы сузить список, уточни запрос "
                "(например, символом гена вместо полного названия белка)."
            )
        st.markdown(
            f"**Запрос неоднозначный:** найдено {len(candidates)} похожих записей."
            + total_note
            + " Выбери запись — карточка соберётся по её коду UniProt."
        )
        labels = {
            c["accession"]: (
                f"{c['accession']} — {c['название_белка']} · "
                f"{c['организм_отображение']}"
            )
            for c in candidates
        }
        choice = st.selectbox(
            "Запись в UniProt",
            ["— выбери запись —"] + list(labels),
            format_func=lambda a: labels.get(a, a),
            key="выбор_записи",
        )
        if choice == "— выбери запись —":
            st.info("Выбери запись из списка, и здесь появится карточка.")
            return
        accession = choice

    # Уровень уверенности определяет СПОСОБ сопоставления — передаём его
    # как есть, а факт ручного выбора идёт отдельным флагом и уровень НЕ
    # подменяет. Иначе почти любой поиск получал бы «ВЫСОКУЮ»: кандидатов
    # больше одного почти всегда (у insulin их 5435), человек ткнул в первую
    # строчку — и слабый текстовый поиск выглядел бы как точное определение.
    # Ручной выбор честно дописывается к причине — см. core.assess_confidence.
    user_confirmed = len(candidates) > 1
    search_way = way_by_accession.get(accession, "")

    with st.spinner("Опрашиваю UniProt, NCBI и RCSB PDB… обычно это 5–15 секунд"):
        card = cached_card(accession, search_way, user_confirmed)

    # ----------------------------------------------------------------------
    # Ничего не нашлось
    # ----------------------------------------------------------------------
    if "uniprot" not in card["источники"]:
        st.error("Ничего не нашлось в UniProt — без этого источника карточку не собрать.")
        for problem in card["проблемы_источников"]:
            st.markdown(f"- {problem}")
        return

    up = card["источники"]["uniprot"]
    conf = card.get("уверенность", {})
    level = conf.get("уровень", "")
    level_key = level.lower()
    reason = conf.get("причина", "")

    # ----------------------------------------------------------------------
    # Шапка карточки: название, организм, уверенность
    # ----------------------------------------------------------------------
    st.divider()

    # Научное название организма — латинское, по правилам биологической
    # номенклатуры пишется курсивом; бытовое русское название — перед ним.
    ru_org = core.ORGANISM_RU.get(up["организм"])
    organism_md = (
        f"{ru_org} (*{up['организм']}*)" if ru_org else f"*{up['организм']}*"
    )

    left, right = st.columns([2, 1])
    with left:
        st.subheader(up["название_белка"])
        st.markdown(
            f"**UniProt:** [{up['accession']}]({up['ссылка_на_запись']}) "
            f"({up['id_записи']})  \n"
            f"**Гены:** {', '.join(g for g in up['гены'] if g) or '—'}  \n"
            f"**Организм:** {organism_md}  \n"
            f"**Статус проверки записи:** {core.entry_status_display(up.get('статус_проверки'))}"
        )
    with right:
        # st.metric предназначен для чисел — текстовый уровень в нём выглядит
        # сломанным. Поэтому цветной бейдж; расшифровка причины — плашкой ниже.
        if level_key.startswith(("высок", "high")):
            badge = ":green[ВЫСОКАЯ]"
        elif level_key.startswith(("средн", "medium")):
            badge = ":orange[СРЕДНЯЯ]"
        else:
            badge = ":red[НИЗКАЯ]"
        st.markdown(f"**Уверенность сопоставления**\n\n### {badge}")

    # Уровни уверенности ядро пишет по-русски ("высокая"/"средняя"/"низкая"),
    # проверяем начало слова, чтобы не зависеть от формы.
    if level_key.startswith(("высок", "high")):
        st.success(f"**Высокая уверенность.** {reason}")
    elif level_key.startswith(("средн", "medium")):
        st.info(f"**Средняя уверенность.** {reason}")
    elif level_key.startswith(("низк", "low")):
        st.warning(f"**Низкая уверенность** — запись может быть не тем, что ты искала. {reason}")

    # ----------------------------------------------------------------------
    # Метрики
    # ----------------------------------------------------------------------
    m1, m2, m3 = st.columns(3)
    m1.metric("Длина белка", f"{up['длина_последовательности_ак']} а.к.")
    m2.metric("3D-структур в PDB", up["число_pdb_структур"])
    m3.metric(
        "Хромосома (NCBI)",
        card["источники"].get("ncbi_gene", {}).get("хромосома") or "—",
    )

    # ----------------------------------------------------------------------
    # Секции источников
    # ----------------------------------------------------------------------
    with st.expander("UniProt — что это за белок", expanded=True):
        st.markdown(f"**Функция:** {linkify_pubmed(up['функция'])}")
        st.markdown(f"**Начало последовательности:** {up['начало_последовательности']}")
        lit = up.get("литература") or []
        if lit:
            with st.expander(f"Литература ({len(lit)} статей в записи UniProt)"):
                for i, pub in enumerate(lit, 1):
                    st.markdown(f"{i}. [PubMed:{pub['pubmed_id']}]({pub['ссылка']})")
                st.caption(
                    "Статьи, которыми UniProt подкрепляет данные записи. "
                    "Полный список — на странице записи UniProt (ссылка в шапке карточки)."
                )

    ng = card["источники"].get("ncbi_gene")
    if ng:
        with st.expander("NCBI Gene — что известно о гене", expanded=True):
            st.markdown(
                f"**Символ гена:** {ng['символ']} (Gene ID [{ng['gene_id']}]({ng['ссылка_на_запись']}))  \n"
                f"**Расположение:** хромосома {ng['хромосома'] or '—'} "
                f"({ng['расположение_в_геноме'] or '—'})  \n"
                f"**Другие названия:** {ng['другие_названия'] or '—'}  \n"
                f"**Описание:** {ng['описание'] or '(нет описания)'}  \n"
                f"**Кратко:** {ng['краткое_описание']}"
            )

    nf = card["источники"].get("ncbi_protein_fallback")
    if nf:
        with st.expander("NCBI Protein — запасной источник (записи гена не нашлось)"):
            st.markdown(
                f"**GenBank ID:** [{nf['genbank_id']}]({nf['ссылка_на_запись']})  \n"
                f"**Заголовок:** {nf['заголовок']}  \n"
                f"**Организм:** {core.organism_display(nf['организм'])}  \n"
                f"**Надёжность:** {nf['надёжность']}"
            )

    pdb = card["источники"].get("pdb")
    if pdb and pdb["структуры"]:
        with st.expander("RCSB PDB — 3D-структуры", expanded=True):
            st.caption(
                f"Известно структур по данным UniProt: {pdb['всего_структур_в_uniprot']}, "
                f"показаны первые {pdb['показано_примеров']}."
            )
            for s in pdb["структуры"]:
                resolution = (
                    f"{s['разрешение_ангстрем']} Å"
                    if s["разрешение_ангстрем"] is not None
                    else "нет (не применимо)"  # у NMR-структур разрешения не бывает
                )
                st.markdown(
                    f"**[{s['pdb_id']}]({s['ссылка']})** — {s['название'] or '(без названия)'}  \n"
                    f"метод: {s['метод'] or '—'} · разрешение: {resolution}"
                )
    elif pdb:
        with st.expander("RCSB PDB — 3D-структуры"):
            st.write("Структуры для этого белка не удалось загрузить.")
    else:
        with st.expander("RCSB PDB — 3D-структуры"):
            st.write("В UniProt для этого белка не указано 3D-структур.")

    # Интерактивная 3D-модель — отдельным блоком под списком структур:
    # вращается мышкой, по умолчанию выключена (не замедляет карточку).
    # Условие — по полному списку ID из UniProt: блок работает, даже если
    # запрос деталей структур в RCSB не удался (файл качается по коду).
    if up.get("pdb_id_список"):
        render_3d_block(pdb, up)

    if card["проблемы_источников"]:
        with st.expander("Что не удалось получить (честно)"):
            for problem in card["проблемы_источников"]:
                st.markdown(f"- {problem}")

    # ----------------------------------------------------------------------
    # Скачивание карточки
    # ----------------------------------------------------------------------
    st.divider()
    # Имя файла — по коду UniProt: в собранной карточке он есть всегда,
    # запасные ветки не нужны (раньше здесь был недостижимый re.sub — только
    # путал при чтении).
    file_stem = up.get("accession") or "запрос"
    c1, c2 = st.columns(2)
    c1.download_button(
        "Скачать карточку (.txt)",
        data=core.render_text(card) + "\n",
        file_name=f"card_{file_stem}.txt",
        mime="text/plain",
        use_container_width=True,
    )
    c2.download_button(
        "Скачать данные (.json)",
        data=json.dumps(card, ensure_ascii=False, indent=2),
        file_name=f"card_{file_stem}.json",
        mime="application/json",
        use_container_width=True,
    )


# ---------------------------------------------------------------------------
# Вкладка 2. Массовая обработка
# ---------------------------------------------------------------------------

def show_batch_results(rows):
    """Сводка + таблица с фильтром + скачивание CSV."""
    summary = core.summarize_batch(rows)

    st.subheader("Сводка")
    cols = st.columns(6)
    cols[0].metric("Всего", summary["всего"])
    cols[1].metric("Найдено", summary["найдено"])
    cols[2].metric("Высокая уверенность", summary["высокая"])
    cols[3].metric("Средняя", summary["средняя"])
    cols[4].metric("Низкая", summary["низкая"])
    not_found = summary["не найдено"] + summary["ошибок"]
    cols[5].metric("Не найдено / ошибки", not_found)

    notes = []
    if summary["неоднозначных"]:
        notes.append(
            f"{summary['неоднозначных']} неоднозначных запросов — для них взята первая "
            "запись по приоритету; проверь колонку «примечание»."
        )
    if summary["низкая"]:
        notes.append(
            f"{summary['низкая']} строк с низкой уверенностью — их стоит проверить вручную "
            "перед использованием в анализе."
        )
    for note in notes:
        st.warning(note)

    df = pd.DataFrame(rows)

    st.subheader("Результаты")
    filter_choice = st.selectbox(
        "Показать строки",
        ["все", "только низкая уверенность", "только не найденные и ошибки"],
        key="batch_filter",
    )
    if filter_choice == "только низкая уверенность":
        df = df[df["уверенность"] == "НИЗКАЯ"]
    elif filter_choice == "только не найденные и ошибки":
        df = df[df["результат"] != "найдено"]

    st.dataframe(df, hide_index=True, use_container_width=True)

    stamp = time.strftime("%Y%m%d_%H%M")
    st.download_button(
        "Скачать CSV",
        data=df.to_csv(index=False).encode("utf-8-sig"),
        file_name=f"biolink_batch_{stamp}.csv",
        mime="text/csv",
        use_container_width=False,
    )
    st.caption(
        "CSV сохранён в кодировке UTF-8 с BOM — открывается в Excel без кракозябр. "
        "Колонка «ссылка_uniprot» — кликабельная ссылка на запись."
    )


def render_batch_tab():
    st.markdown(
        "Вставь список названий (через запятую или с новой строки) или загрузи файл "
        "(в CSV берётся первая колонка). Для каждого названия программа соберёт "
        "карточку из трёх баз и добавит строку в таблицу с уровнем уверенности "
        "сопоставления.\n\n"
        "**Честное предупреждение:** в массовом режиме программа не спрашивает по "
        "каждому названию, как вкладка «Поиск», — если записей несколько, берётся "
        "первая по приоритету. Спорные случаи (совпадение по синониму или обычным "
        "текстовым поиском) помечаются «неоднозначно» и попадают в счётчик сводки. "
        "Когда же программа просто выбрала человека среди организмов с тем же "
        "именем гена (TP53 есть и у крысы, и у мыши), в «примечании» стоит "
        "спокойная пометка — это нормальная работа приоритета, а не риск. "
        "Сомнительные строки проверяются по колонке «примечание»."
    )

    uploaded = st.file_uploader(
        "Загрузить файл (.txt / .csv)",
        type=["txt", "csv"],
    )
    raw_uploaded = ""
    if uploaded is not None:
        raw_bytes = uploaded.getvalue()
        try:
            raw_uploaded = raw_bytes.decode("utf-8-sig")
        except UnicodeDecodeError:
            try:
                # Excel на Windows часто сохраняет CSV в cp1251.
                raw_uploaded = raw_bytes.decode("cp1251")
            except UnicodeDecodeError:
                st.error("Не удалось прочитать файл как UTF-8 или cp1251.")
                raw_uploaded = ""

    raw_text = st.text_area(
        "…или вставь список здесь",
        value="",
        height=150,
        placeholder="insulin, TP53, p53\nGFP, P01308",
        disabled=bool(raw_uploaded),
    )
    if raw_uploaded:
        st.caption(
            "Файл загружен, поэтому ручной ввод отключён. Чтобы снова печатать "
            "в поле, убери файл — крестик на загрузчике выше."
        )

    raw = raw_uploaded or raw_text
    is_csv = uploaded is not None and uploaded.name.lower().endswith(".csv")
    if is_csv:
        terms, duplicates = core.parse_batch_input(raw, csv_mode=True)
    else:
        terms, duplicates = core.parse_batch_input(raw)

    allowed_to_start = bool(terms)
    if not terms:
        if raw:
            st.info("Не нашлось ни одного названия — проверь список.")

    if terms:
        est_seconds = len(terms) * EST_SECONDS_PER_CARD
        key_state = "с ключом NCBI (быстрее)" if core.NCBI_API_KEY else "без ключа NCBI"
        st.info(
            f"Загружено названий: **{len(terms)}**"
            + (f" (повторов пропущено: {duplicates})" if duplicates else "")
            + f". Оценка времени: **{format_duration(est_seconds)}** ({key_state})."
        )

        if len(terms) > 300:
            st.warning(
                f"Список большой: {len(terms)} названий при {EST_SECONDS_PER_CARD} сек "
                f"на карточку — это {format_duration(est_seconds)}. Вкладку закрывать нельзя, "
                "пока обработка не закончится. Турбо-режим для больших списков "
                "(пакетный ID mapping) добавим следующим шагом."
            )
            if not st.checkbox("Подтверждаю: я готова ждать", key="batch_confirm"):
                allowed_to_start = False

    start = False
    if allowed_to_start:
        start = st.button("Начать обработку", type="primary")

    # Пока новый прогон не запущен — показываем результаты прошлого. Без этого
    # любой клик (фильтр, «Скачать CSV») перерисовывал страницу в пустоту:
    # Streamlit при каждом клике перезапускает скрипт целиком, а результаты
    # рисовались только один раз — сразу после конца прогона.
    if not start:
        saved = st.session_state.get("batch_rows")
        if saved:
            st.divider()
            st.caption(
                "Результаты последнего прогона — можно фильтровать и скачивать. "
                "Кнопка «Начать обработку» выше запустит новый прогон."
            )
            show_batch_results(saved)
        return

    cache = st.session_state.setdefault("batch_cache", {})
    # Ограничение размера: кэш нужен, чтобы повторные названия не дёргали
    # базы заново, но расти бесконечно он не должен (несколько прогонов по
    # сотне названий — и вкладка тяжелеет). Храним не более 5000 последних
    # записей: dict помнит порядок вставки, выкидываем самые старые.
    if len(cache) > 5000:
        for old_term in list(cache)[: len(cache) - 4000]:
            cache.pop(old_term, None)
    rows = []
    # Сохраняем список в session_state сразу: он будет пополняться после
    # каждой карточки, и если пользователь случайно кликнет посреди долгого
    # прогона, частичный результат не пропадёт (Streamlit прерывает скрипт
    # при любом взаимодействии и перезапускает его).
    st.session_state["batch_rows"] = rows
    started = time.time()

    progress = st.progress(0.0)
    status = st.empty()
    table_placeholder = st.empty()
    st.caption(
        "Прервать прогон можно в любой момент — просто кликни по любому элементу "
        "страницы (например, поменяй фильтр): обработанные к этому моменту строки "
        "сохранятся, и их можно будет скачать."
    )

    for i, term in enumerate(terms):
        if term in cache:
            row = cache[term]
        else:
            try:
                # детали_pdb=False: в таблицу нужны только ЧИСЛО структур
                # (оно уже приходит из UniProt) и уверенность, а не названия
                # и разрешения структур. Минус до трёх запросов RCSB на белок
                # — список любой длины обрабатывается примерно вдвое быстрее.
                card = core.build_card(term, детали_pdb=False)
            except Exception as exc:  # сеть, таймаут — обработка не должна падать
                card = {
                    "запрос_пользователя": term,
                    "источники": {},
                    "проблемы_источников": [f"ошибка: {exc}"],
                }
            row = core.batch_row_from_card(card)
            if row.get("результат") != "ошибка":
                # Ошибочные строки НЕ кэшируем: «ошибка» из-за моргнувшей
                # сети не должна приклеиваться к названию до конца сессии —
                # повторный прогон обязан повторить запрос по-настоящему.
                cache[term] = row
        rows.append(row)

        done = i + 1
        avg = (time.time() - started) / done
        remaining = avg * (len(terms) - done)
        progress.progress(done / len(terms))
        status.markdown(
            f"**Обработано {done} из {len(terms)}** — осталось примерно "
            f"{format_duration(remaining)}. Текущий запрос: `{term}`"
        )
        if done % 5 == 0 or done == len(terms):
            table_placeholder.dataframe(
                pd.DataFrame(rows), hide_index=True, use_container_width=True
            )
        # Паузы между карточками здесь не нужны: вежливые интервалы уже
        # соблюдаются внутри ядра (time.sleep(ncbi_delay()) перед каждым
        # запросом к NCBI), а между двумя NCBI-запросами всегда вклинивается
        # запрос к UniProt (0.3–1 с). Остаточный риск очень длинных прогонов:
        # если сервер всё же ответит 429, ядро повторит запрос с нарастающей
        # паузой, а совсем не удавшиеся строки честно помечаются «ошибка».

    st.success("Готово! Результаты ниже — можно отфильтровать и скачать в CSV.")

    show_batch_results(rows)


def clear_batch_results():
    if st.session_state.get("batch_rows"):
        st.divider()
        if st.button("Очистить результаты массовой обработки"):
            st.session_state.pop("batch_rows", None)
            st.session_state.pop("batch_cache", None)
            st.rerun()


tab_search, tab_batch = st.tabs(["Поиск карточки", "Массовая обработка"])

# Порядок блоков совпадает с порядком вкладок — так код читается легче.
with tab_search:
    render_search_tab()

with tab_batch:
    render_batch_tab()
    clear_batch_results()

st.divider()
st.caption(
    "Данные: [UniProt](https://www.uniprot.org) · [NCBI](https://www.ncbi.nlm.nih.gov) · "
    "[RCSB PDB](https://www.rcsb.org). Программа не анализирует белок сама — она честно "
    "объединяет данные баз и показывает уверенность сопоставления."
)

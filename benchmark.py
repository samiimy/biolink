# -*- coding: utf-8 -*-
"""Мини-бенчмарк точности программы (prototype_semifinale).

Идея из ревью: чтобы проект стал исследовательской работой, нужно измерить,
насколько точно программа сопоставляет запрос пользователя с записью в UniProt.

Как это честно устроено:
1. GOLDEN — «эталонный набор»: 52 известных белка + их правильные UniProt-коды.
   НЕТ веры на слово: перед тестом каждый эталон проверяется отдельным
   запросом к UniProt (живой ли код, тот ли ген, тот ли организм).
2. Для каждого запроса вызывается НАСТОЯЩАЯ функция поиска из
   prototype_semifinale (fetch_uniprot) и настоящая assess_confidence —
   бенчмарк не подменяет логику.
3. Метрики:
   - Полнота (recall): для скольких запросов вообще найдена запись.
   - Точность (accuracy): точное попадание / попадание в семейство / ошибка.
   - Калибровка уверенности: совпадает ли заявленная уверенность с реальной
     точностью (ВЫСОКАЯ должна ошибаться реже, чем НИЗКАЯ).
"""

import csv
import json
import sys
import time

from prototype_semifinale import SESSION, fetch_uniprot, assess_confidence, HEADERS

# ---------------------------------------------------------------- эталоны
# (запрос, категория, ожидаемый accession, ген для проверки, организм для
#  проверки или None=человек, допустимое семейство генов или None)
G = "человек"
GOLDEN = [
    # --- A. точные UniProt-коды ---
    ("P01308", "A:код",       "P01308", "INS",    G,    None),
    ("P04637", "A:код",       "P04637", "TP53",   G,    None),
    ("P02768", "A:код",       "P02768", "ALB",    G,    None),
    ("P69905", "A:код",       "P69905", "HBA1",   G,    None),
    ("P68871", "A:код",       "P68871", "HBB",    G,    None),
    # --- B. основные имена генов ---
    ("INS",    "B:ген",       "P01308", "INS",    G,    None),
    ("TP53",   "B:ген",       "P04637", "TP53",   G,    None),
    ("BRCA1",  "B:ген",       "P38398", "BRCA1",  G,    None),
    ("BRCA2",  "B:ген",       "P51587", "BRCA2",  G,    None),
    ("HBA1",   "B:ген",       "P69905", "HBA1",   G,    None),
    ("HBB",    "B:ген",       "P68871", "HBB",    G,    None),
    ("ALB",    "B:ген",       "P02768", "ALB",    G,    None),
    ("CFTR",   "B:ген",       "P13569", "CFTR",   G,    None),
    ("APOE",   "B:ген",       "P02649", "APOE",   G,    None),
    ("ESR1",   "B:ген",       "P03372", "ESR1",   G,    None),
    ("EGFR",   "B:ген",       "P00533", "EGFR",   G,    None),
    ("KRAS",   "B:ген",       "P01116", "KRAS",   G,    None),
    ("MYC",    "B:ген",       "P01106", "MYC",    G,    None),
    ("TNF",    "B:ген",       "P01375", "TNF",    G,    None),
    ("IL6",    "B:ген",       "P05231", "IL6",    G,    None),
    ("BDNF",   "B:ген",       "P23560", "BDNF",   G,    None),
    ("DMD",    "B:ген",       "P11532", "DMD",    G,    None),
    ("HTT",    "B:ген",       "P42858", "HTT",    G,    None),
    ("F8",     "B:ген",       "P00451", "F8",     G,    None),
    ("PTPRC",  "B:ген",       "P08575", "PTPRC",  G,    None),
    # --- C. синонимы генов и имён белков ---
    ("p53",        "C:синоним", "P04637", "TP53",    G, None),
    ("HER2",       "C:синоним", "P04626", "ERBB2",   G, None),
    ("CD45",       "C:синоним", "P08575", "PTPRC",   G, None),
    ("OCT4",       "C:синоним", "Q01860", "POU5F1",  G, None),
    ("p21",        "C:синоним", "P38936", "CDKN1A",  G, None),
    ("p16",        "C:синоним", "P42771", "CDKN2A",  G, None),
    ("MDR1",       "C:синоним", "P08183", "ABCB1",   G, None),
    ("Rb",         "C:синоним", "P06400", "RB1",     G, None),
    ("Huntingtin", "C:синоним", "P42858", "HTT",     G, None),
    ("Epo",        "C:синоним", "P01588", "EPO",     G, None),
    # --- D. текстовые имена белков (тут честно ожидаем низкую уверенность) ---
    # Для неоднозначных слов эталон — правильное СЕМЕЙСТВО генов, не один код.
    ("insulin",        "D:текст", "P01308", "INS",   G, None),
    ("hemoglobin",     "D:текст", None, "HB",        G, ("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBE1", "HBM", "HBQ1")),
    ("albumin",        "D:текст", "P02768", "ALB",   G, None),
    ("transferrin",    "D:текст", "P02787", "TF",    G, None),
    ("dystrophin",     "D:текст", "P11532", "DMD",   G, None),
    ("catalase",       "D:текст", "P04040", "CAT",   G, None),
    ("trypsin",        "D:текст", "P07477", "PRSS1", G, None),
    ("pepsin",         "D:текст", "P00790", "PGA",   G, None),
    ("growth hormone", "D:текст", None, "GH",     G, ("GH1", "GH2", "CSH1", "CSH2")),
    ("keratin",        "D:текст", None, "KRT",       G, "KRT"),
    ("collagen",       "D:текст", None, "COL",       G, "COL"),
    ("amylase",        "D:текст", None, "AMY",       G, "AMY"),
    ("myosin",         "D:текст", None, "MY",        G, "MYH|MYO|MYL"),
    # --- E. нечеловеческие белки ---
    ("GFP",       "E:не-человек", "P42212", "GFP",     "Aequorea victoria", None),
    ("lacZ",      "E:не-человек", "P00722", "lacZ",    "Escherichia coli",  None),
    ("BSA",       "E:не-человек", "P02769", "ALB",     "Bos taurus",        None),
    ("ovalbumin", "E:не-человек", "P01012", "OVAL",    "Gallus gallus",     None),
]


def verify_golden():
    """Независимая проверка эталонов: живой ли код, тот ли ген/организм."""
    accs = sorted({r[2] for r in GOLDEN if r[2]})
    ok = {}
    # UniProt позволяет одним запросом получить много записей по accession.
    for i in range(0, len(accs), 25):
        chunk = accs[i:i + 25]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        url = "https://rest.uniprot.org/uniprotkb/search"
        r = SESSION.get(url, params={
            "query": q, "format": "json", "size": 100,
            "fields": "accession,gene_primary,organism_name",
        }, headers=HEADERS, timeout=30)
        r.raise_for_status()
        for e in r.json().get("results", []):
            genes = [g.get("geneName", {}).get("value")
                     for g in (e.get("genes") or []) if g]
            org = ((e.get("organism") or {}).get("scientificName") or "")
            ok[e["primaryAccession"]] = (genes[0] if genes else None, org)
        time.sleep(0.4)

    verified, dropped = [], []
    for row in GOLDEN:
        query, cat, acc, gene, org, fam = row
        if acc is None:  # эталон-семейство: нечего проверять по коду
            verified.append(row)
            continue
        info = ok.get(acc)
        if not info:
            dropped.append((row, "код не найден в UniProt (эталон ошибочен)"))
        elif info[0] != gene:
            dropped.append((row, f"код жив, но ген {info[0]!r} != эталон {gene!r}"))
        elif (org == G and not info[1].startswith("Homo sapiens")):
            dropped.append((row, f"ожидали человека, а организм {info[1]!r}"))
        else:
            verified.append(row)
    return verified, dropped


def judge(found, row):
    """Точное попадание / семейство / ошибка / не найдено."""
    query, cat, acc, gene, org, fam = row
    if found is None:
        return "не найдено"
    f_acc = found.get("accession")
    f_genes = []
    for g in (found.get("гены") or []):
        v = g.get("value") if isinstance(g, dict) else g
        if v:
            f_genes.append(v)
    f_org = (found.get("организм") or "")
    want_human = (org == G)
    is_human = f_org.startswith("Homo sapiens")

    if acc and f_acc == acc and ((not want_human) or is_human):
        return "точно"
    if fam:  # эталон-семейство
        if want_human and not is_human:
            return "ошибка"  # семейство верное, но не тот организм
        for gv in f_genes:
            if isinstance(fam, tuple) and gv in fam:
                return "семейство"
            if isinstance(fam, str):
                prefixes = fam.split("|")
                if any(gv.startswith(p) for p in prefixes):
                    return "семейство"
    return "ошибка"


def main():
    print("=== ШАГ 1: верификация эталонов ===")
    verified, dropped = verify_golden()
    for row, why in dropped:
        print(f"  [ВЫБЫЛ] {row[0]!r} ({row[1]}): {why}")
    print(f"Эталонов проверено: {len(verified)} прошло, {len(dropped)} выбыло\n")

    rows = []
    for n, row in enumerate(verified, 1):
        query, cat = row[0], row[1]
        print(f"[{n}/{len(verified)}] {cat} {query!r} ...")
        try:
            found, err = fetch_uniprot(query)
            if err:
                print(f"  !! поиск не удался: {err}")
        except Exception as e:  # noqa: BLE001 — бенчмарк не должен падать
            print(f"  !! исключение: {e}")
            found = None
        conf = assess_confidence(found)[0] if found else "—"
        verdict = judge(found, row)
        got_acc = (found or {}).get("accession")
        way = (found or {}).get("способ_поиска", "")
        if isinstance(way, dict):
            way = way.get("gene", "") or str(way)
        rows.append({
            "категория": cat, "запрос": query, "найдено": got_acc or "",
            "эталон": row[2] or ("семейство " + str(row[5])),
            "вердикт": verdict, "уверенность": conf, "способ": way,
        })
        print(f"  -> {got_acc} | {verdict} | уверенность: {conf}")
        time.sleep(0.5)

    # ------------------------------------------------------------- метрики
    total = len(rows)
    found_rows = [r for r in rows if r["вердикт"] != "не найдено"]
    exact = [r for r in rows if r["вердикт"] == "точно"]
    family = [r for r in rows if r["вердикт"] == "семейство"]
    wrong = [r for r in rows if r["вердикт"] == "ошибка"]

    def pct(x):
        return f"{100 * x / total:.0f}%" if total else "—"

    print("\n=== ИТОГИ ===")
    print(f"Всего запросов:       {total}")
    print(f"Полнота (найдено):    {len(found_rows)}/{total} = {pct(len(found_rows))}")
    print(f"Точное попадание:     {len(exact)}/{total} = {pct(len(exact))}")
    print(f"Верное семейство:     {len(family)}/{total} = {pct(len(family))}")
    print(f"Ошибки:               {len(wrong)}/{total} = {pct(len(wrong))}")

    print("\n--- По категориям ---")
    cats = sorted({r["категория"] for r in rows})
    by_cat = {}
    for c in cats:
        cr = [r for r in rows if r["категория"] == c]
        good = sum(r["вердикт"] in ("точно", "семейство") for r in cr)
        fnd = sum(r["вердикт"] != "не найдено" for r in cr)
        by_cat[c] = (len(cr), fnd, good)
        print(f"{c:14s}: найдено {fnd}/{len(cr)}, верно {good}/{len(cr)}")

    print("\n--- Калибровка уверенности (главное!) ---")
    calib = {}
    for lvl in ("высокая", "средняя", "низкая"):
        lr = [r for r in rows if r["уверенность"] == lvl]
        if not lr:
            continue
        good = sum(r["вердикт"] in ("точно", "семейство") for r in lr)
        calib[lvl] = (good, len(lr), 100 * good / len(lr))
        print(f"{lvl:8s}: верно {good}/{len(lr)} = {100 * good / len(lr):.0f}%")

    # ------------------------------------------------------------- вывод
    with open("бенчмарк-результаты.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    with open("бенчмарк-итоги.json", "w", encoding="utf-8") as f:
        json.dump({
            "total": total,
            "recall": len(found_rows),
            "exact": len(exact),
            "family": len(family),
            "wrong": len(wrong),
            "по_категориям": {k: v for k, v in by_cat.items()},
            "калибровка": calib,
            "выбывшие_эталоны": [(r[0], w) for r, w in dropped],
        }, f, ensure_ascii=False, indent=2)

    print("\nСохранено: бенчмарк-результаты.csv, бенчмарк-итоги.json")


if __name__ == "__main__":
    main()

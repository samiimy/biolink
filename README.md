# BioLink — one protein card from three scientific databases

*English summary below; the full documentation is in Russian (scroll down).*

BioLink is a small Streamlit web app that takes a gene or protein name (`insulin`, `TP53`, `p53`) or a UniProt accession (`P01308`) and builds a single card from three open databases:

- **UniProt**: protein function, organism, amino-acid sequence;
- **NCBI Gene**: chromosome location, official description, aliases (falls back to NCBI Protein when no Gene record exists, e.g. GFP);
- **RCSB PDB**: solved 3D structures (method, resolution) plus an interactive 3D viewer in the browser (3Dmol.js via py3Dmol).

The card is exported as human-readable text and as JSON, so BioLink works as a thin *linking layer* between databases instead of opening three sites by hand.

## What makes it more than a wrapper: calibrated match confidence

Name-based lookup is ambiguous (`p53`, `HER2` match many records). BioLink never silently picks the first hit in interactive mode: it shows the candidates and lets the user choose. Every card also reports **how** the record was matched:

- exact UniProt accession → high confidence
- query equals the primary gene name → high
- query is a gene synonym → medium
- free-text search → low

A manual choice by the user is recorded as a separate fact, not as "high confidence", so that the scale keeps its meaning.

**Batch mode** processes a list of names (text, `.txt` or `.csv`), picks records by transparent priority rules (similar to `--pick` / `--flag_pick` in Ensembl VEP), flags ambiguous rows, and exports a CSV.

## Evaluation

`benchmark.py` runs the real lookup function against a gold set of 52 well-known proteins (each gold accession is re-verified live against UniProt; 2 were dropped). Results, September 2026:

- recall: 50/50 (100%)
- exact accession match: 38/50 (76%); right gene family: 4/50; wrong: 8/50
- **confidence calibration**: "high" is correct in 28/29 cases (97%), "medium" in 5/7 (71%), "low" in 9/14 (64%)

The error rate therefore increases precisely where the application reports a weak match. The principal result of the project is calibration rather than raw accuracy.

## Quick start

```bash
pip install -r requirements.txt   # streamlit, requests, pandas, py3Dmol
streamlit run app.py
```

An NCBI API key is optional (put it in `.streamlit/secrets.toml`, see the `.example` template); without it the app simply respects NCBI rate limits. A console version: `python3 prototype_semifinale.py insulin`.

## Limitations

This is a prototype, not a predictive AI system: it links existing records and does not infer biology. It does not check that a PDB structure belongs to the same isoform as the UniProt record, does not compare text descriptions across databases, and does not version data. Batch mode relies on priority rules rather than guarantees, which is why ambiguous rows are flagged.

Related heavy-weight tools that solve data integration at scale: BioMart, InterMine, GeneCards. BioLink is a lightweight learning prototype of the same idea, whose distinctive part is the benchmarked confidence scale.

The choice of data sources, the card design, the confidence scale and its benchmark were developed manually; the code was written with the use of AI coding assistants. License: MIT.

---

# BioLink — карточка белка из трёх научных баз по одному запросу

Веб-приложение на Python (Streamlit). На вход подаётся название гена или белка (например, `insulin`, `TP53`, `p53`) либо идентификатор UniProt (например, `P01308`). Приложение запрашивает три открытые научные базы данных:

- **UniProt** — функция белка, организм, аминокислотная последовательность;
- **NCBI Gene** — данные о гене: хромосомная локализация, официальное описание, альтернативные названия;
- **RCSB PDB** — наличие расшифрованных 3D-структур (метод, разрешение, ссылки) и **интерактивная 3D-модель в браузере** (движок 3Dmol.js через py3Dmol).

Полученные данные объединяются в одну карточку, доступную в браузере, в текстовом файле и в формате JSON. Приложение выполняет роль связующего слоя между базами: вместо ручного просмотра трёх сайтов достаточно одного запроса.

## Режимы работы

### 1. Поиск карточки

Один запрос формирует карточку с разделами по каждой базе, метриками и шкалой уверенности сопоставления.

**Неоднозначные запросы.** Если по названию (например, `p53` или `HER2`) в UniProt найдено несколько записей, приложение не выбирает первую из них автоматически: оно показывает список кандидатов (идентификатор, название, организм), и карточка формируется после выбора пользователя. Для точного идентификатора UniProt (`P01308`) карточка формируется сразу. Организм указывается общеупотребительным названием с научным (латинским) названием в скобках: «Человек (*Homo sapiens*)».

**Шкала уверенности** является ключевой частью проекта: карточка указывает, каким способом запись была сопоставлена с запросом.

| Способ сопоставления | Уверенность |
|---|---|
| точный идентификатор UniProt | высокая |
| совпадение с основным названием гена | высокая |
| совпадение с синонимом гена | средняя |
| текстовый поиск | низкая |

Если запись выбрана пользователем из списка, уровень уверенности определяется способом поиска, а в карточке отдельно отмечается, что выбор подтверждён вручную. Ручной выбор учитывается как отдельный факт, а не как уровень шкалы: в противном случае почти любой поиск получал бы высокий уровень, и шкала утратила бы смысл.

### 2. Пакетная обработка

На вход подаётся список названий (через запятую, с новой строки или файлом `.txt`/`.csv`). Приложение последовательно формирует карточки и выводит таблицу: найденный идентификатор, организм, число структур, уровень уверенности, примечание. Результат выгружается в CSV (UTF-8 с BOM, корректно открывается в Excel). В ходе обработки отображается прогресс; обработку можно прервать, обработанные строки сохраняются.

**Ограничение пакетного режима:** в отличие от режима поиска, приложение не запрашивает выбор для каждого названия; при нескольких записях выбирается первая по приоритету. Спорные случаи (совпадение по синониму или текстовому поиску) отмечаются как неоднозначные и учитываются в сводке. Выбор записи человека среди нескольких организмов с тем же названием гена (TP53 есть также у крысы и мыши) отмечается в примечании как штатная работа приоритета. Аналогичный подход применяется в крупных инструментах аннотации (например, `--pick` / `--flag_pick` в Ensembl VEP): автоматический выбор по прозрачным правилам с явной пометкой спорных строк.

## Запуск

```bash
pip install -r requirements.txt      # streamlit, requests, pandas, py3Dmol
streamlit run app.py
```

Приложение откроется в браузере по адресу `http://localhost:8501`.

**Доступ с телефона в той же сети Wi-Fi:** определите IP-адрес компьютера (Windows: `ipconfig`, строка «IPv4-адрес») и откройте `http://<IP>:8501`. При первом запуске Windows запросит разрешение брандмауэра для частных сетей.

**Публикация в интернете:** бесплатный хостинг [Streamlit Community Cloud](https://streamlit.io/cloud) разворачивает приложение непосредственно из GitHub-репозитория.

### Ключ NCBI API (необязательно)

Без ключа NCBI допускает около 3 запросов в секунду, с бесплатным ключом — до 10, что ускоряет пакетную обработку. Получение ключа: ncbi.nlm.nih.gov → регистрация → Settings → API Key Management → Create an API Key.

Ключ указывается в файле `.streamlit/secrets.toml` (файл не входит в репозиторий; шаблон — `.streamlit/secrets.toml.example`):

```toml
[ncbi]
api_key = "ваш_ключ"
```

При отсутствии файла приложение работает медленнее, соблюдая паузы между запросами к NCBI.

### Консольная версия

Ядро (`prototype_semifinale.py`) запускается без браузера (требуются Python 3 и `requests`):

```bash
python3 prototype_semifinale.py insulin
```

Создаются файлы `card_<запрос>.txt` (для чтения) и `card_<запрос>.json` (машиночитаемый формат для передачи в другие программы). Консольная версия работает в автоматическом режиме и выбирает запись с наивысшим приоритетом, что обеспечивает воспроизводимость бенчмарка.

## Оценка точности: бенчмарк

`benchmark.py` — исследовательская часть проекта. Эталонный набор включает 52 известных белка с правильными идентификаторами UniProt; перед тестом каждый эталон проверяется запросом к UniProt. Для каждого запроса вызывается рабочая функция поиска ядра без подмены логики. Измеряются полнота, точность и **калибровка уверенности**, то есть действительно ли записи с высоким уровнем ошибочны реже, чем с низким.

Результаты на актуальной версии ядра (`prototype_semifinale`, сентябрь 2026):

| Метрика | Значение |
|---|---|
| Эталонов | 50 (2 исключены при проверке) |
| Найдено | 50/50 = 100% |
| Точное совпадение | 38/50 = 76% |
| Верное семейство генов | 4/50 = 8% |
| Ошибки | 8/50 = 16% |
| Калибровка уровня «высокая» | 28/29 = **97%** |
| Калибровка уровня «средняя» | 5/7 = 71% |
| Калибровка уровня «низкая» | 9/14 = 64% |

Основная метрика — калибровка: при уровне «высокая» сопоставление верно в 97% случаев (28 из 29 на эталонной выборке), а доля ошибок возрастает именно там, где приложение предупреждает о слабом сопоставлении. Повторный прогон: `python3 benchmark.py` (требуется доступ в интернет, около 2–3 минут).

## Назначение и границы применимости

BioLink **не является** предсказательной ИИ-системой. Приложение не выполняет предсказаний: оно находит и объединяет существующие записи трёх баз по названию гена или белка. Система, предсказывающая поведение генетической конструкции в живой клетке, представляет собой отдельную и значительно более сложную задачу; её обучение требует данных в едином формате, и связующий слой такого типа является первым шагом к их накоплению.

## Известные ограничения

Проект является прототипом. Открытые ограничения:

- **Отсутствует проверка согласованности между базами.** Приложение не проверяет, что найденная в PDB структура относится к той же изоформе белка, что и запись UniProt (идентификаторы структур берутся из ссылок UniProt).
- **Отсутствует семантическая нормализация.** Описание функции из UniProt и описание гена из NCBI выводятся рядом без сопоставления.
- **Отсутствует версионирование данных.** Карточка отражает состояние баз на момент запроса.
- **Пакетный режим основан на приоритетах, а не на гарантии.** При нескольких совпадениях выбирается запись человека с основным названием гена, спорные строки отмечаются. Это документированная практика (см. Ensembl VEP `--pick`), но не гарантия, поэтому предусмотрены колонка «примечание» и шкала уверенности.
- **Темп запросов в пакетном режиме определяется ядром.** Пауза выдерживается перед каждым запросом к NCBI, а между двумя запросами к NCBI всегда выполняется запрос к UniProt. При отказе сервера (код 429) запрос повторяется с увеличивающейся паузой; необработанные строки отмечаются в таблице как ошибочные.
- **Отсутствует ускоренный режим для больших списков.** Запросы выполняются последовательно (около 3 секунд на карточку); пакетное сопоставление идентификаторов для тысяч названий запланировано.

## Технические детали

- **UniProt** возвращает все необходимые поля, включая список структур PDB, одним запросом.
- **NCBI Gene** содержит записи только для организмов с аннотированным геномом. Для инсулина человека (INS) запись подробная, а для GFP (медуза *Aequorea victoria*) записи в NCBI Gene нет; это особенность базы, а не ошибка приложения. В таком случае приложение переключается на NCBI Protein и отмечает это как запасной вариант.
- **RCSB PDB** поддерживает пакетные запросы (GraphQL), однако прототип запрашивает структуры по одной; у распространённых белков структур сотни (у инсулина — 382, у GFP — 654), поэтому показываются первые примеры.
- **3D-модель** отображается движком 3Dmol.js (тем же, что используется в stmol); файл структуры загружается с `files.rcsb.org` и кэшируется.
- **Ключ NCBI API** необязателен: без него пауза между запросами к NCBI составляет около 0.4 с, с ключом — около 0.15 с. Ключ хранится в `.streamlit/secrets.toml`, который не входит в репозиторий (шаблон прилагается).
- У части структур PDB (например, полученных методом ЯМР) отсутствует поле «разрешение», так как оно определено только для рентгеноструктурных данных; в карточке это указывается явно.

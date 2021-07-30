# v2,v3 (STAR BES)

Данный код позволяет измерить эллиптический и треугольный потоки в столкновениях ядер золота при энергиях BES STAR, используя метод плоскости события.

## Содержание:

I. [Устновка](#Устновка) \
II. [подготовка данных](#ПодготовкаДанных) \
III. [Обработка событий](#EventProcessing) \

## I. Устновка <a name="Устновка"></a>

На кластере NICA

```bash
cd /scratch2/${USER}
mkdir STAR
cd STAR
git clone https://github.com/DemanovAE/BES.git
cd BES
```

Не забудьте добавить библиотеки ROOT в свою среду, используя thisroot.sh
В терминале кластера NICA:

```sh
source /opt/fairsoft/bmn/may18p1/bin/thisroot.sh
```

Компилирование читалки данных:
(собранная читалка для данных femtoDst лежит по пути `/scratch2/demanov/STAR/BES/StFemtoEvent/libStFemtoDst.so`)

```bash
cd StFemtoEvent
make
```
Измените пути в `set_env.sh`: изменить `ST_FEMTO_DST_INC_DIR` на стандартный. Это каталог, в котором хранится `libStFemtoDst.so`.

Установка проекта с помощью CMake. (находясь в директории ./BES)

```bash
mkdir build
cd build
cmake ../macro
make
```

## II. Подготовка данных <a name="ПодготовкаДанных"></a>

Используйте `./scripts/GenerateLists.sh` для создания списков файлов:

```bash
. GenerateLists.sh FEMTODST_DIR N_FILES_IN_LIST
```
где FEMTODST_DIR - путь к каталогу с femtoDst.root файлами. И N_FILES_IN_LIST обозначает максимальное количество femtoDst.root файлов в каждом листе. Базовый пример: `. GenerateLists.sh /scratch2/parfenov/StData/27gev/run1/ 100`

Полученные списки файлов будут в BES/lists/


## III. Обработка событий <a name="EventProcessing"></a>

Интерактивный режим:
```bash
./FemtoDstAnalyzer_PID -i inFile -o outFile -m WorkMode -g Energy
```
1. `inFile` - root файл или лист с файлами
2. `outFile` - выходной файл. Для измерения потоков нужно провести 3 прогонки данных. Выходной файл после 1 прогонки `NoRe_27GeV.root`, после второй - `Re_27GeV.root` и конечный файл `Flow_27GeV.root`.
3. `WorkMode` - указывает стадию анализа.

        | WorkMode        | Описание |
        | --------------- | ----------- |
        | QA              | в разработке
        | raw             | Первая прогонка данных. На данном этапе набираются данные для корекции Q-векторов, а именно дял процедуры реценренинга
        | rec             | Вторая прогонка данных. На данном этапе набираются данные для корекции угла плоскости события, а именно дял процедуры флатенинга
        | flow            | Третья прогонка данных. На данном этапе измеряются v2, v3, и расрешение.
4. `Energy` - энергия анализируемых данных.

Пример запуска:
```bash
./FemtoDstAnalyzer_PID -i ../lists/lists27GeV/StRun15.list -o ./NoRe_27GeV.root -m raw -g 27
```

Для отправки задач на класстер NICA используются следующий bash скрипт - `/scripts/start_pid_nica.sh INPUT_FILELIST_DIR INPUT_WORKMODE INPUT_ENERGY`
```sh
cd /scripts
. start_pid_nica.sh /scratch2/$USER/STAR/BES/lists/lists27GeV raw 27
```
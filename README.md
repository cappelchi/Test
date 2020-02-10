# Test
Решение тестового задания.

Atlas_git.ipynb - полное решение с визуализацией действий

Atlas_.py - скрипт решения

ALL.merged.eu.22chr.4pop.5.Q - Q файл

Минимальное оборудование:
8 core cpu,
160 RAM,
400 HDD

Имеем 404 образца и 78076746 SNPs, каждый SNPs занимает 4 байта в оперативной памяти, следовательно расчетный инстанс:
404 х 78076746 х 4 / (1024 ^ 3) = 117.5 ГБ оперативной памяти под размещение и ещё ~ 50% от это под инференс, итого ~ 175 ГБ.

![Analyse](https://github.com/cappelchi/Test/blob/master/img/comp_viz.png)

![компоненты](https://github.com/cappelchi/Test/blob/master/img/component.png)

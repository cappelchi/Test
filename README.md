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
404 х 78076746 х 4 / (1024 ^ 3) = 117.5 ГБ оперативной памяти под размещение и ещё + 50% под инференс, итого ~ 175 ГБ.

UPD: 5-мерая PCA проекция европейских популяций с помощью UMAP на плоскость:
![UMAP](https://raw.githubusercontent.com/cappelchi/Test/master/img/umap.gif)

[Интерактивный график](https://plot.ly/~cappelchi/13/)

<img src="https://raw.githubusercontent.com/cappelchi/Test/master/img/umap1.png" width="800" height="600" />

Получается 2 четких финских кластера:

<img src="https://raw.githubusercontent.com/cappelchi/Test/master/img/umap5.png" width="800" height="600" />

Есть небольшой GBR кластер плюс несколько локальных TSi кластеров.

<img src="https://raw.githubusercontent.com/cappelchi/Test/master/img/umap11.png" width="800" height="600" />


Остальные кластеры смешанные:
<img src="https://raw.githubusercontent.com/cappelchi/Test/master/img/umap19.png" width="800" height="600" />

<img src="https://raw.githubusercontent.com/cappelchi/Test/master/img/umap30.png" width="800" height="600" />

<img src="https://raw.githubusercontent.com/cappelchi/Test/master/img/umap39.png" width="800" height="600" />

![Analyse](https://github.com/cappelchi/Test/blob/master/img/comp_viz._g.png)
[Play](https://chart-studio.plot.ly/~cappelchi/10)

![c1](https://github.com/cappelchi/Test/blob/master/img/components.png)

![c2](https://github.com/cappelchi/Test/blob/master/img/components2.png)

![c3](https://github.com/cappelchi/Test/blob/master/img/components3.png)

![c4](https://github.com/cappelchi/Test/blob/master/img/components4.png)

![c5](https://github.com/cappelchi/Test/blob/master/img/components5.png)

![c6](https://github.com/cappelchi/Test/blob/master/img/components6.png)

![c7](https://github.com/cappelchi/Test/blob/master/img/components7.png)

![c8](https://github.com/cappelchi/Test/blob/master/img/components8.png)

![c9](https://github.com/cappelchi/Test/blob/master/img/components9.png)

![c10](https://github.com/cappelchi/Test/blob/master/img/components10.png)

![c11](https://github.com/cappelchi/Test/blob/master/img/components11.png)

![c12](https://github.com/cappelchi/Test/blob/master/img/components12.png)

![c13](https://github.com/cappelchi/Test/blob/master/img/components13.png)

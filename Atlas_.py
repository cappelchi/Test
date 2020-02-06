
# coding: utf-8

# Общий план задания
# 1. Скачать [http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) файлы соответствующие хромосомам 1-22.
# 2. Перевести с помощью plink https://www.cog-genomics.org/plink/1.9/ из vcf в bed\bim\fam
# 3. Оставить в данных только SNP с помощью plink.
# 4. Слить полученные файлы в один набор bed\bim\fam
# 5. Изменить fam файл, чтобы в нем была информация о популяциях. Соответствие образец-популяция можно найти на сайт 1000 геномов.
# 6. Сделать новый набор данных, чтобы в нем были только европейские популяции.
# 7. http://software.genetics.ucla.edu/admixture/ - с помощью этой программы сделать анализ по 5 компонентам (базовый запуск можно найти в документации).
# 
# Как результат ожидается следующее:
# 
# а. Текстовый файл, какие команды использовались для каждого шага (чтобы можно было воспроизвести).
# 
# б.Q файл по итогам 7 шага
# 
# в. Графики разложения отельных индивидов на компоненты

# #### 1. Скачать http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ файлы соответствующие хромосомам 1-22.
# в google colaboratory скачал:
# 
# !echo http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{1..22}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi | xargs -P22 -n1 wget -nv



!wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{str(i)}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz


# #### 2. Перевести с помощью plink https://www.cog-genomics.org/plink/1.9/ из vcf в bed\bim\fam
# #### 3. Оставить в данных только SNP с помощью plink.

# Устанавливаем PLINK
!wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip'
!unzip plink_linux_x86_64_20200121.zip -d plink/

for chr_num in range(1, 23):
    #Распаковываем хромосому в bed-bim-fam файлы
!plink/plink --vcf ALL.chr{str(chr_num)}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed --out ALL.chr{str(chr_num)}
    #Оставляем только SNP
 !plink/plink --bfile ALL.chr{str(chr_num)} --snps-only --make-bed --out ALL.chr{str(chr_num)}
    #Удаляем дубликаты одинаковых SNPs
!cut -f2 ALL.chr{str(chr_num)}.bim | sort | uniq -d > dup{str(chr_num)}.snps
!echo найдено дупликатов в chr{str(chr_num)} хромосоме: 
!wc -l dup{str(chr_num)}.snps
!echo ' - удаляю...'
!plink/plink --bfile ALL.chr{str(chr_num)} --exclude dup{str(chr_num)}.snps --make-bed --out ALL.chr{str(chr_num)}

# #### 4. Слить полученные файлы в один набор bed\bim\fam

#Создаем список всех файлов для объединения в один набор
#!rm merge_list.txt
for i in range(2, 23):
    with open('merge_list.txt', 'a') as mlist:
        mlist.write(f'ALL.chr{str(i)}.bed ALL.chr{str(i)}.bim ALL.chr{str(i)}.fam \n')


!cat merge_list.txt

#Пробуем объединить все файлы в один набор
!plink/plink --threads 6 --bfile ALL.chr1 --merge-list merge_list.txt --make-bed --out ALL.merged


# Найдены трёхаллельные варианты
# Нужно их удалить из файлов
#Удаляем трёхаллельные варианты
for chr_num in range(1, 23):
    !plink/plink --bfile ALL.chr{str(chr_num)} --exclude ALL.merged-merge.missnp --make-bed --out ALL.chr{str(chr_num)}

#Заново сливаем все файлы в один набор
!plink/plink --threads 6 --memory 24576 --bfile ALL.chr1 --merge-list merge_list.txt --make-bed --out ALL.merged


# #### 5. Изменить fam файл, чтобы в нем была информация о популяциях. Соответствие образец-популяция можно найти на сайт 1000 геномов.

#Скачиваем таблицу, где можно найти соответствие образец-популяция.
!wget 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx'

#Устанавливаем и импортируем необходимые библиотеки для работы с таблицей
!pip install pandas
!pip install xlrd

import pandas as pd

#Читаем лист "Sample info"
df = pd.read_excel('20130606_sample_info.xlsx', sheet_name = 'Sample Info', header = 0, usecols = ['Sample', 'Population', 'Population Description'])

#Таблица соответствий
for pop, desc in zip(df['Population'].unique(), df['Population Description'].unique()):
    print (f'{pop} - {desc}')

#Создаем словарь с числовой кодировкой популяций
#Европейские кодируем числом больше > 100
pop_cod = {}
for cnt, samp in enumerate(df.Population.unique()):
    if samp in ['GBR', 'FIN', 'IBS', 'TSI']:
        pop_cod[samp] = cnt + 3 + 100
    else:
        pop_cod[samp] = cnt + 3

#Кодирую популяцию каждого образца в число по ID
df['Popnum'] = 0
ind_cod = {}
for cnt in range(len(df)):
    df.iloc[cnt, 3] = pop_cod[df.Population[cnt]]
    ind_cod[df.iloc[cnt,0]] = pop_cod[df.Population[cnt]]

#Считываем fam - файл
df_fam = pd.read_csv('ALL.merged.fam', sep = ' ', header = None, index_col = False)

#Меняем обозначение популяции на число из словаря
for cnt in range(len(df_fam)):
    df_fam.iloc[cnt, 5] = ind_cod[df_fam.iloc[cnt,0]]

#Перезаписываем fam - файл
df_fam.to_csv('ALL.merged.fam', sep = ' ', header = False, index = False, )

# #### 6. Сделать новый набор данных, чтобы в нем были только европейские популяции.

# Создаём список европейских популяций

df_keep = df_fam.loc[df_fam[5] > 100][1]

print ('Количество европейских образцов = ', len(df_keep))

df_keep.to_csv('keep.fam', sep = ' ', header = False, index = False, )


# Делаем набор из европейских популяций
!plink/plink --bfile ALL.merged --keep-fam keep.fam --make-bed --out ALL.merged.eu


# #### 7. http://software.genetics.ucla.edu/admixture/ - с помощью этой программы сделать анализ по 5 компонентам (базовый запуск можно найти в документации).
!wget http://software.genetics.ucla.edu/admixture/binaries/admixture_linux-1.3.0.tar.gz
!tar -xvf admixture_linux-1.3.0.tar.gz


# Запускаем на 4 итерации.
!adm/admixture -C 4 ALL.merged.eu.bed 5 -j6


# ##### Q файл по итогам 7 шага

# In[41]:


get_ipython().system('head ALL.merged.eu.5.Q')


# ###### Графики разложения отельных индивидов на компоненты
df_Q = pd.read_csv('ALL.merged.eu.5.Q', names = ['comp1', 'comp2', 'comp3', 'comp4', 'comp5'], sep = ' ', header = None, index_col = False)
df_efam = df_fam = pd.read_csv('ALL.merged.europ.fam', sep = ' ', header = None, index_col = False)
df_Q.set_index(df_efam[0])

!pip install plotly.express

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

cod_pop = dict(zip(pop_cod.values(), pop_cod.keys()))

sub_names = []
for cnt in range(len(df_efam)):
    sub_names.append(f'{df_efam.iloc[cnt,0]} - {cod_pop[df_efam.iloc[cnt, 5]]}')

start = 0
end = 40
cols = 4
fig = make_subplots(rows = 10, 
                    cols = cols,
                    subplot_titles = sub_names[start:end],                    
                    #vertical_spacing = 0.08,
                    #horizontal_spacing = None
                   )
for i in range(40):
    fig.append_trace(go.Bar(
        #x = ['comp1', 'comp2', 'comp3', 'comp4', 'comp5'],
        y = df_Q.iloc[i, :],
        hovertext = df_Q.iloc[i, :]
    ), row = i // cols + 1, col = i % cols + 1)
    #print(i // 4 + 1, i % 4 + 1, df_Q.iloc[i, :])
#fig.update_layout(height = 300 + len(df_Q)//4 + 1, width=1000, title_text="Stacked subplots")
fig.update_layout(height = 100 * 10, 
                  width=200 * cols, 
                  title_text="Графики разложения отельных индивидов на компоненты",
                  showlegend = False
                 )
fig.show()


# In[57]:


start = 40
end = 80
cols = 4
fig = make_subplots(rows = 10, 
                    cols = cols,
                    subplot_titles = sub_names[start:end],                    
                    #vertical_spacing = 0.08,
                    #horizontal_spacing = None
                   )
for i in range(40):
    fig.append_trace(go.Bar(
        #x = ['comp1', 'comp2', 'comp3', 'comp4', 'comp5'],
        y = df_Q.iloc[i, :],
        hovertext = df_Q.iloc[i, :]
    ), row = i // cols + 1, col = i % cols + 1)
    #print(i // 4 + 1, i % 4 + 1, df_Q.iloc[i, :])
#fig.update_layout(height = 300 + len(df_Q)//4 + 1, width=1000, title_text="Stacked subplots")
fig.update_layout(height = 100 * 10, 
                  width=200 * cols, 
                  title_text="Графики разложения отельных индивидов на компоненты",
                  showlegend = False
                 )
fig.show()


# In[58]:


start = 80
end = 120
cols = 4
fig = make_subplots(rows = 10, 
                    cols = cols,
                    subplot_titles = sub_names[start:end],                    
                    #vertical_spacing = 0.08,
                    #horizontal_spacing = None
                   )
for i in range(40):
    fig.append_trace(go.Bar(
        #x = ['comp1', 'comp2', 'comp3', 'comp4', 'comp5'],
        y = df_Q.iloc[i, :],
        hovertext = df_Q.iloc[i, :]
    ), row = i // cols + 1, col = i % cols + 1)
    #print(i // 4 + 1, i % 4 + 1, df_Q.iloc[i, :])
#fig.update_layout(height = 300 + len(df_Q)//4 + 1, width=1000, title_text="Stacked subplots")
fig.update_layout(height = 100 * 10, 
                  width=200 * cols, 
                  title_text="Графики разложения отельных индивидов на компоненты",
                  showlegend = False
                 )
fig.show()


# In[59]:


start = 120
end = 160
cols = 4
fig = make_subplots(rows = 10, 
                    cols = cols,
                    subplot_titles = sub_names[start:end],                    
                    #vertical_spacing = 0.08,
                    #horizontal_spacing = None
                   )
for i in range(40):
    fig.append_trace(go.Bar(
        #x = ['comp1', 'comp2', 'comp3', 'comp4', 'comp5'],
        y = df_Q.iloc[i, :],
        hovertext = df_Q.iloc[i, :]
    ), row = i // cols + 1, col = i % cols + 1)
    #print(i // 4 + 1, i % 4 + 1, df_Q.iloc[i, :])
#fig.update_layout(height = 300 + len(df_Q)//4 + 1, width=1000, title_text="Stacked subplots")
fig.update_layout(height = 100 * 10, 
                  width=200 * cols, 
                  title_text="Графики разложения отельных индивидов на компоненты",
                  showlegend = False
                 )
fig.show()


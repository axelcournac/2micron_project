import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.proportion import proportions_chisquare
from scipy.stats import chi2_contingency
import pandas as pd


# fichier CSV dans un DataFrame
file_path = "/home/axel/Bureau/2micron_plasmid_PROJECT/stabilities_2024/stability_cinetic2.csv"  # Remplacez par le chemin vers votre fichier CSV
# file_path = "/home/axel/Bureau/2micron_plasmid_PROJECT/stabilities_2024/stability_cinetic1and2_merged.csv"  
df = pd.read_csv(file_path)

# Renommer les colonnes en fonction des points temporels
df.columns = ['souche', 'replicat', 'condition', 'time_0h', 'time_6h', 'time_12h', 'time_24h']

# DataFrame en dictionnaire avec les listes pour chaque colonne
data = {
    'souche': df['souche'].tolist(),
    'replicat': df['replicat'].tolist(),
    'condition': df['condition'].tolist(),
    'time_0h': df['time_0h'].tolist(),
    'time_6h': df['time_6h'].tolist(),
    'time_12h': df['time_12h'].tolist(),
    'time_24h': df['time_24h'].tolist()
}


df = pd.DataFrame(data)

# Remplacement des noms des souches
df['souche'] = df['souche'].replace({1218: '5toA', 1068: 'deltaREP1', 1069: 'deltaSTB', 1070: 'WT'})

# Calcul des ratios G418/YPD pour chaque réplicat, souche et point temporel
ratios = pd.DataFrame()
for time in ['time_0h', 'time_6h', 'time_12h', 'time_24h']:
    # Sélectionner les données pour YPD et G418 séparément
    ypd_data = df[df['condition'] == 'YPD'][['souche', 'replicat', time]].rename(columns={time: 'YPD'})
    g418_data = df[df['condition'] == 'G418'][['souche', 'replicat', time]].rename(columns={time: 'G418'})
    
    # Fusionner les deux DataFrames pour calculer les ratios
    merged_data = pd.merge(ypd_data, g418_data, on=['souche', 'replicat'])
    merged_data['ratio'] = merged_data['G418'] / merged_data['YPD']
    merged_data['time'] = time
    ratios = pd.concat([ratios, merged_data])

# Réordonner les souches dans l'ordre désiré (WT, 5toA, deltaREP1, deltaSTB)
ratios['souche'] = pd.Categorical(ratios['souche'], categories=['WT', '5toA', 'deltaREP1', 'deltaSTB'], ordered=True)

# if we want to remove time point 0 :
ratios = ratios[ratios['time'] != "time_0h"]    

# Box or bar plot des ratios pour chaque temps et souche
plt.figure(figsize=(12, 8))
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

plt.rc('font', **font)

# sns.boxplot(x='time', y='ratio', hue='souche', data=ratios, showfliers=False)

sns.barplot(x='time', y='ratio', hue='souche',  
            data=ratios,errorbar="sd", capsize=0.05, dodge=True, linewidth=1.5, 
            edgecolor='black')

ratios = ratios.reset_index(drop=True)

sns.stripplot(x='time', y='ratio', hue='souche', data=ratios,
    jitter=True, linewidth=0.5, dodge=True, legend=False)

# plt.title("Ratios G418/YPD par souche et temps")
plt.ylabel("Proportion of cells having the plasmid")
plt.xlabel("")

# Calcul des p-values (test du chi-deux sur les proportions) entre les souches WT et 5toA pour chaque point temporel
time_points = ['time_0h', 'time_6h', 'time_12h', 'time_24h']
for time in time_points:
    # Sélection des données pour WT et 5toA
    g418_WT = df[(df['souche'] == 'WT') & (df['condition'] == 'G418')][time].sum()
    ypd_WT = df[(df['souche'] == 'WT') & (df['condition'] == 'YPD')][time].sum()
    g418_5toA = df[(df['souche'] == '5toA') & (df['condition'] == 'G418')][time].sum()
    ypd_5toA = df[(df['souche'] == '5toA') & (df['condition'] == 'YPD')][time].sum()
    
    # Test du chi-deux sur les proportions uniquement sur YPD
    count = [g418_WT, g418_5toA]
    total = [ypd_WT, ypd_5toA]  # Le total est uniquement basé sur YPD

    # Test du chi-deux sur les proportions
    # first function 
    chi2_stat, p_value1, _ = proportions_chisquare(count , total)    
    print(f'p = {p_value1:.2e}')
    
    res= chi2_contingency([count , total])   # allow to have more precise p-values
    chi2_stat, p_value = res.statistic, res.pvalue
    print(f'p = {p_value:.2e}')
    
    # Affichage des p-values sur le plot au-dessus de chaque point temporel
    plt.text(time_points.index(time), 0.6, f'p = {p_value:.2e}', ha='center', fontsize=12)

# Ajustement des limites de l'axe Y pour placer correctement les p-values
plt.ylim(0, ratios['ratio'].max() + 0.2)
plt.show()














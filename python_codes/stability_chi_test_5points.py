import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.proportion import proportions_chisquare
from scipy.stats import chi2_contingency
import pandas as pd


# Représentation du nouveau jeu de données (avec 5 réplicats pour WT et 5toA)
data = {
    'souche': ['WT', 'WT', 'WT', 'WT', 'WT', '5toA', '5toA', '5toA', '5toA', '5toA'],
    'YPD': [382, 319, 1235, 1180, 577, 640, 1025, 343, 109, 227],
    'G418': [141, 93, 656, 568, 190, 183, 265, 107, 44, 62]
}

df = pd.DataFrame(data)

# Calcul des ratios G418/YPD pour chaque souche
df['ratio'] = df['G418'] / df['YPD']

# Réordonner les souches dans l'ordre désiré (WT, 5toA)
df['souche'] = pd.Categorical(df['souche'], categories=['WT', '5toA'], ordered=True)

# Box plot des ratios pour les deux souches
plt.figure(figsize=(8, 6))
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

# sns.boxplot(x='souche', y='ratio', hue='souche',  data=df, showfliers=False)

sns.barplot(x='souche', y='ratio', hue='souche',  
            data=df,errorbar="sd", capsize=0.1, dodge=False, linewidth=1.5, 
            edgecolor='black')

plt.xlabel("Strain")
# plt.bar(x='souche', height='ratio', data=df, yerr=np.std())

# ratio = ratio.reset_index(drop=True)
sns.stripplot(x='souche', y='ratio', hue='souche', data=df,
    jitter=True, linewidth=0.5, dodge=False, s=10.0)

# plt.title("Ratios G418/YPD par souche")
plt.ylabel("Proportion of cells having the plasmid")

# Calcul des p-values (test du chi-deux sur les proportions) entre les souches WT et 5toA
g418_WT = df[df['souche'] == 'WT']['G418'].sum()
ypd_WT = df[df['souche'] == 'WT']['YPD'].sum()
g418_5toA = df[df['souche'] == '5toA']['G418'].sum()
ypd_5toA = df[df['souche'] == '5toA']['YPD'].sum()

# Test du chi-deux sur les proportions
count = [g418_WT, g418_5toA]
total = [ypd_WT, ypd_5toA]

# Test du chi-deux sur les proportions
# chi2_stat, p_value, _ = proportions_chisquare(count, total)

# other function
res= chi2_contingency([count , total])   # allow to have more precise p-values
chi2_stat, p_value = res.statistic, res.pvalue

# Affichage des résultats de p-value sur le plot
plt.text(0.5, df['ratio'].max() + 0.05, f'p = {p_value:.2e}', ha='center', fontsize=12)

# Ajustement des limites de l'axe Y pour placer correctement la p-value
plt.ylim(0, df['ratio'].max() + 0.1)
plt.legend()
plt.show()



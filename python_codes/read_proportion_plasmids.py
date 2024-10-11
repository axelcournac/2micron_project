import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# Données des mutants
data_5toA = [15.145571768427, 14.9747704008977, 16.9130990926959]
data_WT = [13.7843504560051, 11.5387799142154, 12.8782839847123]

# Moyenne et écart-type pour chaque mutant
mean_5toA = np.mean(data_5toA)
mean_WT = np.mean(data_WT)
std_5toA = np.std(data_5toA)
std_WT = np.std(data_WT)

# Calcul du t-test
t_stat, p_value = stats.ttest_ind(data_WT, data_5toA, equal_var=False)

# Création du bar plot avec barres d'erreur
labels = ['WT', '5toA']
means = [mean_WT, mean_5toA]
errors = [std_WT, std_5toA]
colors = ['#ff7f0e', '#1f77b4']  # WT en orange et 5toA en bleu

plt.figure(figsize=(6, 6))
plt.bar(labels, means, yerr=errors, color=colors, edgecolor='black', linewidth=1.5, capsize=5)

# Ajouter des labels et un titre
plt.xlabel('Mutants', fontsize=12)
plt.ylabel('Mean Ratio Avg Number / Percentage Cells with Plasmid', fontsize=12)
plt.title('Comparison of WT and 5toA Mutants', fontsize=14)

# Afficher la p-value sur le graphique
plt.text(0.5, max(means) + max(errors) * 0.1, f'p-value: {p_value:.4f}', ha='center', fontsize=12, color='red')

# Affichage du graphique
plt.tight_layout()
plt.show()

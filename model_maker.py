# Import necessary libraries
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
from sklearn.svm import SVC
import pandas as pd
import numpy as np
import pickle

# Define the paths to the TSV files
path = "protein_info_completed.tsv"

# Define a dictionary that maps each subcellular location string to a corresponding word from the subcellular_locations list
location_mapping = {
    'mitochondrion': 'mitochondria',
    'mitochondria': 'mitochondria',
    'reticulum': 'endoplasmic reticulum',
    'endoplasmic': 'endoplasmic reticulum',
    'apparatus': 'golgi',
    'golgi': 'golgi',
    'lysosome': 'lysosome',
    'peroxisome': 'peroxisome',
    'chloroplast': 'chloroplast',
    'Mitochondrion': 'mitochondria',
    'Mitochondria': 'mitochondria',
    'Reticulum': 'endoplasmic reticulum',
    'Endoplasmic': 'endoplasmic reticulum',
    'Apparatus': 'golgi',
    'Golgi': 'golgi',
    'Lysosome': 'lysosome',
    'Peroxisome': 'peroxisome',
    'Chloroplast': 'chloroplast',
    'Cytoskeleton': 'cytosolskeleton',
    'cytoskeleton': 'cytosolskeleton',
}


# Create a list of all TSV file paths
tsv_files = [path]
modified_files=[]
# Iterate over each TSV file
for tsv_file in tsv_files:
    # Read in the TSV file as a DataFrame
    df = pd.read_csv(tsv_file, sep='\t')
        # Drop any rows with missing values
    df = df.dropna()
    # Replace values in the "Subcellular location [CC]" column with a single word from the location_mapping dictionary
    for i, row in df.iterrows():
        subcellular_location = row["Subcellular location [CC]"]
        for location_string, location_word in location_mapping.items():
            if location_string in subcellular_location:
                df.at[i, "Subcellular location [CC]"] = location_word
                break  # Stop searching for a match once a location is found
            else:
                df.at[i, "Subcellular location [CC]"] =subcellular_location
    # Save the modified DataFrame as a new TSV file with the same name as the original TSV file, but with "_mod" appended to the end
    new_tsv_file = tsv_file[:-4] + "_mod.tsv"
    df.to_csv(new_tsv_file, sep='\t', index=False)
    modified_files.append(tsv_file[:-4] + "_mod.tsv")
    loca=tsv_file[:-4] + "_mod.tsv"
    break

# Define the paths to the TSV files
path = "protein_info_completed_mod.tsv"


location_mapping = {
    'mitochondrion': 'mitochondria',
    'mitochondria': 'mitochondria',
    'reticulum': 'endoplasmic reticulum',
    'endoplasmic': 'endoplasmic reticulum',
    'apparatus': 'golgi',
    'golgi': 'golgi',
    'lysosome': 'lysosome',
    'peroxisome': 'peroxisome',
    'cytoplasm': 'cytosol',
    'cytoskeleton': 'cytosol',
    'plasma': 'cell membrane',
    'chloroplast': 'chloroplast',
    'Mitochondrion': 'mitochondria',
    'Mitochondria': 'mitochondria',
    'Reticulum': 'endoplasmic reticulum',
    'Endoplasmic': 'endoplasmic reticulum',
    'Apparatus': 'golgi',
    'Golgi': 'golgi',
    'Lysosome': 'lysosome',
    'Peroxisome': 'peroxisome',
    'Cytoplasm': 'cytosol',
    'Cytoskeleton': 'cytosolskeleton',
    'Plasma': 'cell membrane',
    'Chloroplast': 'chloroplast',
    'membrane':'cell membrane',
    'nucleus': 'nucleus',
    'Nucleus': 'nucleus',
    'No location': '',
}

# Create a list of all TSV file paths
tsv_files = [path]
modified_files=[]
# Iterate over each TSV file
for tsv_file in tsv_files:
    # Read in the TSV file as a DataFrame
    df = pd.read_csv(tsv_file, sep='\t')
        # Drop any rows with missing values
    df = df.dropna()
    # Replace values in the "Subcellular location [CC]" column with a single word from the location_mapping dictionary
    for i, row in df.iterrows():
        subcellular_location = row["Subcellular location [CC]"]
        for location_string, location_word in location_mapping.items():
            if location_string in subcellular_location:
                df.at[i, "Subcellular location [CC]"] = location_word
                break  # Stop searching for a match once a location is found
            else:
                df.at[i, "Subcellular location [CC]"] ='No location'
    # Save the modified DataFrame as a new TSV file with the same name as the original TSV file, but with "_mod" appended to the end
    new_tsv_file = tsv_file[:-4] + "_mod.tsv"
    df.to_csv(new_tsv_file, sep='\t', index=False)
    modified_files.append(tsv_file[:-4] + "_mod.tsv")

    
# Read in all the modified TSV files and stack them together into a single DataFrame
#modified_files = [endo_path[:-4] + "_mod.tsv", cell_path[:-4] + "_mod.tsv", mito_path[:-4] + "_mod.tsv", lyso_path[:-4] + "_mod.tsv", vacu_path[:-4] + "_mod.tsv", nucl_path[:-4] + "_mod.tsv"]
dfs = []
for file in modified_files:
    df = pd.read_csv(file, sep='\t')
    dfs.append(df)
stacked_data = pd.concat(dfs, axis=0, ignore_index=True)

# Count how many duplicates
duplicates = stacked_data['Entry'].duplicated()
num_duplicates = duplicates.sum()

# Remove duplicates
data_no_duplicates = stacked_data.drop_duplicates(subset='Entry', keep='first')
duplicates = data_no_duplicates['Entry'].duplicated()
num_duplicates = duplicates.sum()
data = data_no_duplicates
data.to_csv('number_of_occurances.tsv', sep='\t', index=False)


data = data.drop(data[data['Subcellular location [CC]'] == 'No location'].index)

list=['cytosol','mitochondria','endoplasmic reticulum','golgi','chloroplast','lysosome','peroxisome','nucleus','cell membrane']

for i in list:
    val=data["Subcellular location [CC]"].value_counts()[i]
    print(i)
    print(val)

  
dfs = []
for file in modified_files:
    df = pd.read_csv(file, sep='\t')
    dfs.append(df)
stacked_data = pd.concat(dfs, axis=0, ignore_index=True)

# Count how many duplicates
duplicates = stacked_data['Entry'].duplicated()
num_duplicates = duplicates.sum()

print(f'Number of duplicate entries: {num_duplicates}')

# Remove duplicates
data_no_duplicates = stacked_data.drop_duplicates(subset='Entry', keep='first')
duplicates = data_no_duplicates['Entry'].duplicated()
num_duplicates = duplicates.sum()

data = data_no_duplicates
data = data.drop(data[data['Subcellular location [CC]'] == 'No location'].index)

# Convert Sequences to k-mers
def kmer_generator(Sequence, k=3):
    return [Sequence[i:i+k] for i in range(len(Sequence)-k+1)]

data["kmer_Sequence"] = data["Sequence"].apply(kmer_generator)
df=data

# Define the list of words to filter by
words_to_filter = ['chloroplast','lysosome','peroxisome']

# Filter rows to remove any row where any word in the "Subcellular location [CC]" column occurs more than 1000 times
for word in words_to_filter:
    mask = df['Subcellular location [CC]'].str.contains(word)
    count = mask.sum()
    if count > 1000:
        drop_indices = df[mask].index[3000:]
        df = df.drop(drop_indices)
    elif count > 0:
        pass
# Define the list of words to filter by
words_to_filter = ['mitochondria','cell membrane']

# Filter rows to remove any row where any word in the "Subcellular location [CC]" column occurs more than 1000 times
for word in words_to_filter:
    mask = df['Subcellular location [CC]'].str.contains(word)
    count = mask.sum()
    if count > 3000:
        drop_indices = df[mask].index[3500:]
        df = df.drop(drop_indices)
    elif count > 0:
        pass 
# Define the list of words to filter by
words_to_filter = ['cytosol','nucleus']

# Filter rows to remove any row where any word in the "Subcellular location [CC]" column occurs more than 1000 times
for word in words_to_filter:
    mask = df['Subcellular location [CC]'].str.contains(word)
    count = mask.sum()
    if count > 3000:
        drop_indices = df[mask].index[10000:]
        df = df.drop(drop_indices)
    elif count > 0:
        pass 

# Define the list of words to filter by
words_to_filter = ['endoplasmic reticulum','golgi']

# Filter rows to remove any row where any word in the "Subcellular location [CC]" column occurs more than 1000 times
for word in words_to_filter:
    mask = df['Subcellular location [CC]'].str.contains(word)
    count = mask.sum()
    if count > 3000:
        drop_indices = df[mask].index[2850:]
        df = df.drop(drop_indices)
    elif count > 0:
        pass 
    
list=['cytosol','mitochondria','endoplasmic reticulum','golgi','chloroplast','lysosome','peroxisome','nucleus','cell membrane']

for i in list:
    val=df["Subcellular location [CC]"].value_counts()[i]
    print(i)
    print(val)

import csv

list=['cytosol','mitochondria','endoplasmic reticulum','golgi','chloroplast','lysosome','peroxisome','nucleus','cell membrane']

with open('subcellular_location.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Subcellular Location', 'Count'])
    
    for i in list:
        val = df['Subcellular location [CC]'].value_counts()[i]
        writer.writerow([i, val])
        
data=df

# Train a classifier
# Define features and labels
X = data["kmer_Sequence"].str.join(" ")
y = data["Subcellular location [CC]"]


# # Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=45)

# Vectorize the k-mer Sequences
vectorizer = TfidfVectorizer()
X_train_tfidf = vectorizer.fit_transform(X_train)
X_test_tfidf = vectorizer.transform(X_test)

# Encode the labels
label_encoder = LabelEncoder()
y_train_encoded = label_encoder.fit_transform(y_train)
y_test_encoded = label_encoder.transform(y_test)
###clear
# Train the classifier
classifier = MLPClassifier(hidden_layer_sizes=(100,), max_iter=10000, random_state=45)
classifier.fit(X_train_tfidf, y_train_encoded)

# Test the classifier
y_pred = classifier.predict(X_test_tfidf)

# Get the unique labels in the test set
unique_labels = np.unique(np.concatenate((y_test_encoded, y_pred)))

# Output results
# print(classification_report(y_test_encoded, y_pred, labels=unique_labels, target_names=label_encoder.inverse_transform(unique_labels)))
# Output results
print(classification_report(y_test_encoded, y_pred, labels=unique_labels, target_names=label_encoder.inverse_transform(unique_labels), zero_division=1))
info_map=classification_report(y_test_encoded, y_pred, labels=unique_labels, target_names=label_encoder.inverse_transform(unique_labels), zero_division=1)
# Open a text file in write mode
with open("my_file.txt", "w") as file:
    # Write the value of the variable to the file
    file.write(info_map)
# Visualize the results
cm = confusion_matrix(y_test_encoded, y_pred, labels=unique_labels)

plt.imshow(cm, cmap=plt.cm.Blues)
plt.xlabel("Predicted labels")
plt.ylabel("True labels")
plt.xticks(np.arange(len(unique_labels)), label_encoder.inverse_transform(unique_labels), rotation=90)
plt.yticks(np.arange(len(unique_labels)), label_encoder.inverse_transform(unique_labels))
plt.title("Confusion matrix")
plt.colorbar()
plt.show()

# # Generate the classification report
report = classification_report(y_test_encoded, y_pred, labels=unique_labels, target_names=label_encoder.inverse_transform(unique_labels), zero_division=1)

# Convert the report to a DataFrame
report_df = pd.DataFrame.from_dict(classification_report(y_test_encoded, y_pred, labels=unique_labels, target_names=label_encoder.inverse_transform(unique_labels), zero_division=1, output_dict=True))

# Save the DataFrame to a CSV file
report_df.to_csv('classification_report.csv', index=False)
# Save the trained model to a file
with open('trained_model.pkl', 'wb') as f:
    pickle.dump(classifier, f)
    

with open('label_encoder.pkl', 'wb') as f:
    pickle.dump(label_encoder, f)

with open('vectorizer.pkl', 'wb') as f:
    pickle.dump(vectorizer, f)
# print(unique_labels)
np.save('my_array.npy', cm)


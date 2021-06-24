#### script to assess motifs and compare survival of TARGET-ALL study #### 

import pandas as pd
import numpy as np
import scipy.stats as stats
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt

kmf = KaplanMeierFitter()
results = []


##### IMMUNE RECOVERIES data ##### 
path_recoveries = r"/path_to_immune_reads.csv"
df = pd.read_csv(path_recoveries)


###### SURVIVAL DATA #######
path_survival = r"H:/Research/CURRENT-PROJECTS/ALL/all_phase2_target_2018_pub_clinical_data.tsv"
data_surv = pd.read_csv(path_survival, sep="\t")
data_surv = data_surv[["Patient ID",'Overall Survival Days','Overall Survival (Months)', 'Overall Survival Status',]]
data_survival = data_surv.rename(columns = {"Patient ID":"PATIENT_ID",
                                            "Overall Survival Status":"OS_STATUS_CATEGORICAL",
                                            "Overall Survival (Months)":"OS_MONTHS"})
data_survival.loc[data_survival["OS_STATUS_CATEGORICAL"] == "DECEASED", "OS_STATUS"] = 1
data_survival.loc[data_survival["OS_STATUS_CATEGORICAL"] == "LIVING", "OS_STATUS"] = 0


###### clinical data #######
survival = r"/all_phase2_target_2018_pub_clinical_data.tsv"
data_surv = pd.read_csv(survival, sep="\t")
data_clinical = data_surv[["Patient ID","Cell of tumor origin"]]


def translate_cdr3_to_sequence(seq):
    positive = ["K","R"]
    negative = ["D","E"]
    aromatic = ["W","F","Y"]
    hydrophobic = ["A","V","L","I","M"]
    other = ["G","H","C","S","T","N","Q","P","X"]

    def join_seq(seq):
        sep = ""
        converted_seq =  sep.join(seq)
        return converted_seq

    new_seq = []  
    
    for aa in seq:
        if aa in positive:
            new_seq.append("P")
        elif aa in negative:
            new_seq.append("N")
        elif aa in aromatic:
            new_seq.append("A")
        elif aa in hydrophobic:
            new_seq.append("H")
        elif aa in other:
            new_seq.append("X")
        else:
            print("unable to convert {} in {}".format(aa,seq))
     
    converted_seq = join_seq(new_seq)
    
    return converted_seq


def getKmers(sequence, size=5):
    ls_kmers = [sequence[x:x+size] for x in range(len(sequence) - size + 1)]
    return ls_kmers


def find_unique_kmers(df):
    unique = []
    for ls in df.KMERS:
        for kmer in ls:
            if kmer not in unique:
                unique.append(kmer)
    return unique


def filter_df(df,filetype,sampletype,receptor):
    df = df[df["File_type"] == filetype]
    df = df[df["Sample Type"] == sampletype]
    df = df[df["Receptor"] == receptor]
    return df


def generate_MOTIF_group(motif, df):
    patients = []
    df_motif = df[df["Translated_CDR3"].str.contains(motif)]
    if not df_motif.empty:
        patients = df_motif["Case ID"].drop_duplicates().tolist()
    return patients
    

####### GENERATE REMAINING CASE IDS  GROUP ######
def generate_ar_group(filetype,sampletype, patients):
    list_AR = []
    path_ar = (r"H:/Research/CURRENT-PROJECTS/ALL/ALL_all_patients_per_groups.xlsx")
    df_ar = pd.read_excel(path_ar, sheet_name=filetype)
    list_ALL = df_ar["AR - {}".format(sampletype)].dropna().values.tolist() 
    
    list_AR = [patient_ALL for patient_ALL in list_ALL if patient_ALL not in patients]
    return list_AR

def get_repetitions_groups(df,cdr3_count,filetype,sampletype,receptor):
    
    pivot = pd.pivot_table(df, index = ["File_type","Sample Type","Receptor",'Case ID',"CDR3"], aggfunc ='size')
    pivot_reps = pivot[pivot >= cdr3_count]
    
    try:
        patients_with_reps = pivot_reps.loc[filetype,sampletype,receptor].index.get_level_values(0).drop_duplicates().tolist()
        unique_patients_with_reps = []
        unique_patients_with_reps = [patient for patient in patients_with_reps if patient not in unique_patients_with_reps]
    
    except:
        unique_patients_with_reps = []
        #print("No patient recoveries with repetitions for: {} - {} - {}".format(filetype,sampletype,receptor))
    
    return unique_patients_with_reps

def generate_AR_reps(all_patients,motif_patients):
    AR_reps_group = [patient for patient in all_patients if patient not in motif_patients]
    return AR_reps_group


def get_better_survivor(mean_os_group1,mean_os_group2):
    if mean_os_group1 > mean_os_group2:
        better_survivor = "Top50%"
    elif mean_os_group1 < mean_os_group2:
        better_survivor = "Bottom50%"
    elif mean_os_group1 == mean_os_group2:
        better_survivor = "Equal Survival"
    else:
        better_survivor = ""
    return better_survivor

def get_logrank(df1,df2):
    results = logrank_test(df1["OS_MONTHS"], df2["OS_MONTHS"],
                       event_observed_A=df1["OS_STATUS"], event_observed_B=df2["OS_STATUS"])
    return results.p_value


###### GENERATE RESULTS OF ANALYSIS IN LIST FORMAT ##########    
def get_data(filetype,sampletype,BCR,receptor,cdr3_count,kmer,df1,df2):
    
    N1 = len(df1["PATIENT_ID"].to_list())
    N2 = len(df2["PATIENT_ID"].to_list())
    
    if N1 and N2 >= 2:
        
        Median_OS_group1 = df1["OS_MONTHS"].median()
        Median_OS_group2 = df2["OS_MONTHS"].median()
        Mean_OS_group1 = df1["OS_MONTHS"].mean()
        Mean_OS_group2 = df2["OS_MONTHS"].mean()
        pvalue = get_logrank(df1,df2)
        
        higher_survival = get_better_survivor(Mean_OS_group1,Mean_OS_group2)
        
        data = [filetype,sampletype,BCR,receptor,cdr3_count,
                kmer,N1, N2, Median_OS_group1, Median_OS_group2,
                Mean_OS_group1,Mean_OS_group2,pvalue,higher_survival]
        
    else:
        data = [filetype,sampletype,BCR,receptor,cdr3_count,kmer,N1, N2,"", "","","","",""]
        
    return data



#### set kmers length, ie. 5 in this case ####
KMERS = 5
filetype = "wxs"
sampletypes = df["Sample Type"].drop_duplicates().tolist()
receptors = df["Receptor"].drop_duplicates().tolist()

# translate cdr3 to motif
df['Translated_CDR3'] = df.apply(lambda x: translate_cdr3_to_sequence(x['CDR3']), axis=1)

# generate kmers for all receptors
df['KMERS'] = df.apply(lambda x: getKmers(x['Translated_CDR3'], KMERS), axis=1)
unique_kmers = find_unique_kmers(df)

print("# of unique KMERs:",len(unique_kmers))



###### GENERATE KM CURVES ##### 
def generate_km(df1, df2, group1, group2):
    fig, ax = plt.subplots()
    fig = kmf.fit(df1["OS_MONTHS"],
                  event_observed=df1["OS_STATUS"].astype(int), 
                  label= group1)
    ax = kmf.plot()
    fig = kmf.fit(df2["OS_MONTHS"],
                  event_observed=df2["OS_STATUS"].astype(int), 
                  label=group2)
    ax = kmf.plot(ax=ax)
    median_os_group1 = df1["OS_MONTHS"].median()
    median_os_group2 = df2["OS_MONTHS"].median()
    results = logrank_test(df1["OS_MONTHS"], df2["OS_MONTHS"],
                           event_observed_A=df1["OS_STATUS"], event_observed_B=df2["OS_STATUS"])
    print("{} - Median OS: {} months".format(group1,median_os_group1))
    print("{} - Median OS: {} months".format(group2,median_os_group2))
    print("N1:",len(df1.PATIENT_ID))
    print("N2:",len(df2.PATIENT_ID))
    print("p-value:",results.p_value)
    return plt.show()




for sampletype in sampletypes: 
    for receptor in receptors:
        for kmer in unique_kmers:
            
            # filter df based on receptor/sample

            df_filt = filter_df(df,filetype,sampletype,receptor)
            motif_group = generate_MOTIF_group(kmer, df_filt)
            AR_group = generate_ar_group(filetype, sampletype, motif_group)

            df_survival_motif = data_survival[data_survival.PATIENT_ID.isin(motif_group)].dropna().drop_duplicates()
            df_survival_AR = data_survival[data_survival.PATIENT_ID.isin(AR_group)].dropna().drop_duplicates()
            
            pvalue = get_logrank(df_survival_motif, df_survival_AR)
            if pvalue < 0.05 and len(motif_group) > 5:
                print("{} - {} - {} - Motif:{}".format(filetype,sampletype,receptor,kmer))
                print("N - motif:",len(motif_group))
                print("N - AR:",len(AR_group))
                print("p-value:{:4f}".format(pvalue))
                print("")





B_cell_patients = data_clinical.loc[data_clinical["Cell of tumor origin"] == "B Cell ALL"]["Patient ID"].drop_duplicates().tolist()
T_cell_patients = data_clinical.loc[data_clinical["Cell of tumor origin"] == "T Cell ALL"]["Patient ID"].drop_duplicates().tolist()
B_precusor_patients =  data_clinical.loc[data_clinical["Cell of tumor origin"] == "B-Precursor"]["Patient ID"].drop_duplicates().tolist()





def generate_unique_kmers(kmer_len,df):
    KMERS = kmer_len
    df['Translated_CDR3'] = df.apply(lambda x: translate_cdr3_to_sequence(x['CDR3']), axis=1)
    # generate kmers for all receptors
    df['KMERS'] = df.apply(lambda x: getKmers(x['Translated_CDR3'], KMERS), axis=1)
    unique_kmers = find_unique_kmers(df)
    
    return unique_kmers



cols_results = ["Filetype","Sampletype","BCR with +3 copies","Receptor",
                "CDR3 copy count","Motif","N - Motif","N - AR",
                "Median OS - Motif", "Median OS - AR",
                "Mean OS - Motif", "Mean OS - AR",
                "P-value","Better Survivor"]

df_results = pd.DataFrame(results)
df_results.columns = cols_results 



#SCENIC+ PIPELINE from current Seurat+Singac processed datasets (using un-paired data)

export LD_LIBRARY_PATH="/opt/R/3.5.2/lib64/R/lib:$LD_LIBRARY_PATH"
import os
os.environ['R_HOME'] = "/opt/R/3.5.2/lib64/R/"

#Quality control - required for later run SCENIC+
import pybiomart as pbm
dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')
annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].to_numpy(dtype = str)
filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
annot = annot[~filter]
annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
annot = annot[annot.Transcript_type == 'protein_coding']

#############
#R code to generate consensus.bed file for continue this pipeline
df <- data.frame(seqnames=seqnames(combined.peaks),
  starts=start(combined.peaks)-1,
  ends=end(combined.peaks),
  names=c(rep(".", length(combined.peaks))),
  scores=c(rep(".", length(combined.peaks))),
  strands=strand(combined.peaks))

write.table(df, file="combined_peaks_neurogenesis.bed", quote=F, sep="\t", row.names=F, col.names=F)
#############

#Generate fragments dictionary from each of the scATAC_samples
fragments_dict = {
    'basal': '/datos_2/output_basal_atac/_job/outs/fragments.tsv.gz',
    'control_day7': "/datos_2/output_control_7_days_atac/_job/outs/fragments.tsv.gz",
    'control_day13': "/datos_2/output_control_13_days_atac/outs/fragments.tsv.gz",
    'control_day20': "/datos_2/output_control_20_days_atac/_job/outs/fragments.tsv.gz",
    'ba_day7': "/datos_2/output_ba_7_days_atac/_job/outs/fragments.tsv.gz",
    'ba_day13': "/datos_2/output_ba_13_days_atac/outs/fragments.tsv.gz",
    'ba_day20': "/datos_2/output_ba_20_days_atac/_job/outs/fragments.tsv.gz",
    }

#FOR IBEX KAUST SERVER PATHS
fragments_dict = {
    'basal': '/ibex/scratch/projects/c2169/Navarra/Neuro/fragments/fragments_basal.tsv.gz',
    'control_day7': "/ibex/scratch/projects/c2169/Navarra/Neuro/fragments/fragments_c7.tsv.gz",
    'control_day13': "/ibex/scratch/projects/c2169/Navarra/Neuro/fragments/fragments_c13.tsv.gz",
    'control_day20': "/ibex/scratch/projects/c2169/Navarra/Neuro/fragments/fragments_c20.tsv.gz",
    'ba_day7': "/ibex/scratch/projects/c2169/Navarra/Neuro/fragments/fragments_ba7.tsv.gz",
    'ba_day13': "/ibex/scratch/projects/c2169/Navarra/Neuro/fragments/fragments_ba13.tsv.gz",
    'ba_day20': "/ibex/scratch/projects/c2169/Navarra/Neuro/fragments/fragments_ba20.tsv.gz",
    }

#COMPUTE the QC required for laters
from pycisTopic.qc import *
#path_to_regions = '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes/combined_peaks_neurogenesis.bed'

path_to_regions = {'basal': '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/combined_peaks_neurogenesis.bed',
                   'control_day7': '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/combined_peaks_neurogenesis.bed',
                   'control_day13': '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/combined_peaks_neurogenesis.bed',
                   'control_day20': '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/combined_peaks_neurogenesis.bed',
                   'ba_day7': '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/combined_peaks_neurogenesis.bed',
                   'ba_day13': '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/combined_peaks_neurogenesis.bed',
                   'ba_day20': '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/combined_peaks_neurogenesis.bed',
                   }

#FOR IBEX KAUST SERVER PATHS
path_to_regions = {'basal': '/ibex/scratch/projects/c2169/Navarra/Neuro/peaks/peaks_basal.bed',
                   'control_day7': '/ibex/scratch/projects/c2169/Navarra/Neuro/peaks/peaks_c7.bed',
                   'control_day13': '/ibex/scratch/projects/c2169/Navarra/Neuro/peaks/peaks_c13.bed',
                   'control_day20': '/ibex/scratch/projects/c2169/Navarra/Neuro/peaks/peaks_c20.bed',
                   'ba_day7': '/ibex/scratch/projects/c2169/Navarra/Neuro/peaks/peaks_ba7.bed',
                   'ba_day13': '/ibex/scratch/projects/c2169/Navarra/Neuro/peaks/peaks_ba13.bed',
                   'ba_day20': '/ibex/scratch/projects/c2169/Navarra/Neuro/peaks/peaks_ba20.bed',
                   }

metadata_bc, profile_data_dict = compute_qc_stats(
                fragments_dict = fragments_dict,
                tss_annotation = annot,
                stats=['barcode_rank_plot', 'duplicate_rate', 'insert_size_distribution', 'profile_tss', 'frip'],
                label_list = None,
                path_to_regions = path_to_regions,
                n_cpu = 1,
                valid_bc = None,
                n_frag = 100,
                n_bc = None,
                tss_flank_window = 1000,
                tss_window = 50,
                tss_minimum_signal_window = 100,
                tss_rolling_window = 10,
                remove_duplicates = True,
                _temp_dir = '/ibex/scratch/projects/c2169/temp2')
                #_temp_dir = '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/temp/')

os.chdir("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/")
work_dir = "/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/"

#FOR IBEX KAUST SERVER PATHS
os.chdir("/ibex/scratch/projects/c2169/Navarra/Neuro/scenic_outputs")
work_dir = "/ibex/scratch/projects/c2169/Navarra/Neuro/scenic_outputs"

#Number of cells per condition and timepoint
summary = []

for key, df in metadata_bc.items():
    if "_" in key:
        condition, timepoint = key.split("_")
    else:
        condition = key
        timepoint = None
    summary.append((condition, timepoint, df.shape[0]))

import pandas as pd
summary_df = pd.DataFrame(summary, columns=["condition", "timepoint", "n_cells"])
print(summary_df)


sys.stderr = sys.__stderr__  # unsilence stderr

if not os.path.exists(os.path.join(work_dir, 'quality_control')):
    os.makedirs(os.path.join(work_dir, 'quality_control'))

pickle.dump(metadata_bc,
            open(os.path.join(work_dir, 'quality_control/metadata_bc.pkl'), 'wb'))

pickle.dump(profile_data_dict,
            open(os.path.join(work_dir, 'quality_control/profile_data_dict.pkl'), 'wb'))

# Return figure to plot together with other metrics. Figure will be saved as pdf.
qc_filters = {
    'basal': {
        'Log_unique_nr_frag':   [3.8, None],
        'FRIP':                 [0.5, None],
        'TSS_enrichment':       [5, None],
        'Dupl_rate':            [None, None]
    },
    'control_day7': {
        'Log_unique_nr_frag':   [3.8, None],
        'FRIP':                 [0.5, None],
        'TSS_enrichment':       [5, None],
        'Dupl_rate':            [None, None]
    },
    'control_day13': {
        'Log_unique_nr_frag':   [3.8, None],
        'FRIP':                 [0.5, None],
        'TSS_enrichment':       [5, None],
        'Dupl_rate':            [None, None]
    },
    'control_day20': {
        'Log_unique_nr_frag':   [3.8, None],
        'FRIP':                 [0.5, None],
        'TSS_enrichment':       [5, None],
        'Dupl_rate':            [None, None]
    },
    'ba_day7': {
        'Log_unique_nr_frag':   [3.8, None],
        'FRIP':                 [0.5, None],
        'TSS_enrichment':       [5, None],
        'Dupl_rate':            [None, None]
    },
    'ba_day13': {
        'Log_unique_nr_frag':   [3.8, None],
        'FRIP':                 [0.5, None],
        'TSS_enrichment':       [5, None],
        'Dupl_rate':            [None, None]
    },
    'ba_day20': {
        'Log_unique_nr_frag':   [3.8, None],
        'FRIP':                 [0.5, None],
        'TSS_enrichment':       [5, None],
        'Dupl_rate':            [None, None]
    }
}

FRIP_NR_FRAG_filterDict = {}
TSS_NR_FRAG_filterDict = {}
FRIP_NR_FRAG_figDict = {}
TSS_NR_FRAG_figDict = {}
DR_NR_FRAG_figDict={}

from pycisTopic.qc import *
for runID in metadata_bc:
    FRIP_NR_FRAG_fig, FRIP_NR_FRAG_filter=plot_barcode_metrics(metadata_bc[runID],
                                           var_x='Log_unique_nr_frag',
                                           var_y='FRIP',
                                           min_x=qc_filters[runID]['Log_unique_nr_frag'][0],
                                           max_x=qc_filters[runID]['Log_unique_nr_frag'][1],
                                           min_y=qc_filters[runID]['FRIP'][0],
                                           max_y=qc_filters[runID]['FRIP'][1],
                                           return_cells=True,
                                           return_fig=True,
                                           plot=False)
    # Return figure to plot together with other metrics, and cells passing filters
    TSS_NR_FRAG_fig, TSS_NR_FRAG_filter=plot_barcode_metrics(metadata_bc[runID],
                                          var_x='Log_unique_nr_frag',
                                          var_y='TSS_enrichment',
                                          min_x=qc_filters[runID]['Log_unique_nr_frag'][0],
                                          max_x=qc_filters[runID]['Log_unique_nr_frag'][1],
                                          min_y=qc_filters[runID]['TSS_enrichment'][0],
                                          max_y=qc_filters[runID]['TSS_enrichment'][1],
                                          return_cells=True,
                                          return_fig=True,
                                          plot=False)
    # Return figure to plot together with other metrics, but not returning cells (no filter applied for the duplication rate  per barcode)
    DR_NR_FRAG_fig=plot_barcode_metrics(metadata_bc[runID],
                                          var_x='Log_unique_nr_frag',
                                          var_y='Dupl_rate',
                                          min_x=qc_filters[runID]['Log_unique_nr_frag'][0],
                                          max_x=qc_filters[runID]['Log_unique_nr_frag'][1],
                                          min_y=qc_filters[runID]['Dupl_rate'][0],
                                          max_y=qc_filters[runID]['Dupl_rate'][1],
                                          return_cells=False,
                                          return_fig=True,
                                          plot=False,
                                          plot_as_hexbin = True)
    # Barcodes passing filters
    FRIP_NR_FRAG_filterDict[runID] = FRIP_NR_FRAG_filter
    TSS_NR_FRAG_filterDict[runID] = TSS_NR_FRAG_filter
    # Figs
    FRIP_NR_FRAG_figDict[runID] = FRIP_NR_FRAG_fig
    TSS_NR_FRAG_figDict[runID] = TSS_NR_FRAG_fig
    DR_NR_FRAG_figDict[runID]=DR_NR_FRAG_fig

#PLOT QC for eah sample in pdf
runID = 'basal'

print("filter for sample: {runID}")
full_name = work_dir + 'QC_sample' + runID + '.pdf'
fig=plt.figure(figsize=(20, 10), dpi=800)
plt.subplot(1, 3, 1)
img = fig2img(FRIP_NR_FRAG_figDict[runID]) #To convert figures to plot together, see .utils.py
plt.imshow(img)
plt.axis('off')
plt.subplot(1, 3, 2)
img = fig2img(TSS_NR_FRAG_figDict[runID])
plt.imshow(img)
plt.axis('off')
plt.subplot(1, 3, 3)
img = fig2img(DR_NR_FRAG_figDict[runID])
plt.imshow(img)
plt.axis('off')
plt.show()
plt.savefig(full_name)
plt.close()


#############
#Load list of barcodes from R (Seurat to python pipeline)
#bc_passing_filters = pd.read_csv("barcodes_from_filtered_R.csv")

#Al final usar los filtrados de esta manera
bc_passing_filters = dict()
for runID in metadata_bc:
    bc_passing_filters[runID] = list((set(FRIP_NR_FRAG_filterDict[runID]) & set(TSS_NR_FRAG_filterDict[runID])))
    print(f"{len(bc_passing_filters[runID])} barcodes passed filters for sample {runID}")


#Creating a cisTopic object and topic modeling
import pickle
path_to_blacklist= '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/hg38-blacklist.v2.bed'

#FOR IBEX KAUST SERVER PATHS
path_to_blacklist='/ibex/scratch/projects/c2169/Navarra/Neuro/blacklist.v2.bed'

#Next we will create a cisTopic object for each sample, and merge them into a single object
from pycisTopic.cistopic_class import *
cistopic_obj_list=[create_cistopic_object_from_fragments(path_to_fragments=fragments_dict[key],
                                               path_to_regions=path_to_regions[key],
                                               path_to_blacklist=path_to_blacklist,
                                               metrics=metadata_bc[key],
                                               valid_bc=bc_passing_filters[key],
                                               n_cpu=1,
                                               project=key) for key in fragments_dict.keys()]

cistopic_obj = merge(cistopic_obj_list)
print(cistopic_obj)

pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'quality_control/cistopic_obj.pkl'), 'wb'))

#############
#R code 
metadata = integrated@meta.data
barcodes_metadata = sub('.*_', '', rownames(metadata))

metadata["barcodes"] = barcodes_metadata
metadata["batch"] = sub('.*atac_', '', metadata$batch)
write.table(metadata, file="metadata_R.csv", sep= ",", row.names=T, col.names=T)
#############

#############
#Adding cell metadata from Seurat
cell_data = pd.read_csv("metadata_R.csv")
#cistopic makes use of the sample_id to match the correct cell barcodes to the metadata, let's add the sample_id as a suffix to the cell barcodes
cell_data['barcode'] = cell_data['barcodes'] +'___'+ cell_data['batch']
print(cell_data['barcode'][0:5])
cell_data = cell_data.set_index('barcode')
cistopic_obj.add_cell_data(cell_data) #[['MMline']]) cell not contained will be filled with nan

#############
#Detect and remove doublets using Scrublet
import scrublet as scr
scrub = scr.Scrublet(cistopic_obj.fragment_matrix.T, expected_doublet_rate=0.1)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.call_doublets(threshold=0.3)

histogram_scrub = scrub.plot_histogram()
plt.savefig('histogram_scrub.pdf')  

scrublet = pd.DataFrame([scrub.doublet_scores_obs_, scrub.predicted_doublets_], columns=cistopic_obj.cell_names, index=['Doublet_scores_fragments', 'Predicted_doublets_fragments']).T
cistopic_obj.add_cell_data(scrublet)

singlets = cistopic_obj.cell_data[cistopic_obj.cell_data.Predicted_doublets_fragments == False].index.tolist()
cistopic_obj = cistopic_obj.subset(singlets, copy=True)

#######################
#######################
#FROM NOW ON WE USE ONLY THE CELLS FOR THE TRANSITIONS GROUPS OF CELLS A-B-C-D-E

#############
#Subset 
#Remove unwanted clusters / keep wanted clusters of cell based on Seurat annotation from scRNA-seq
our_transition_groups_cells = cistopic_obj.cell_data[cistopic_obj.cell_data['group_A'].isin(['A','B','C','D','E'])].index.tolist()
# Subset cisTopic object
cistopic_obj = cistopic_obj.subset(our_transition_groups_cells, copy=True, split_pattern='-')

#Subset only CONTROL CELLS
#our_transition_groups_cells = cistopic_obj_our_groups_only.cell_data[cistopic_obj_our_groups_only.cell_data['condition'].isin(['control'])].index.tolist()
# Subset cisTopic object
#cistopic_obj_our_groups_only = cistopic_obj_our_groups_only.subset(our_transition_groups_cells, copy=True, split_pattern='-')

#lets try using all cells and use only contrasts we are interested in
cistopic_obj_our_groups_only = cistopic_obj

#######################
#WE ARE WORKING ONLY WITH CONTROL CELLS FROM NEUROGENESIS BRUNCH NOW !!
#######################

#############
#Run topic modeling
import pickle
from pycisTopic.cistopic_class import *
tmp_dir = '/datos_2/temp/'
sys.stderr = open(os.devnull, "w")  # silence stderr
models=run_cgs_models(cistopic_obj_our_groups_only,
                    n_topics=[2,5,10,15,30,45],
                    n_cpu=20,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=None,
                    _temp_dir = os.path.join(tmp_dir + 'ray_spill'))
sys.stderr = sys.__stderr__  # unsilence stderr
#ray.init()
#ray.shutdown()

#############
#Model selection
numTopics = 30
model = evaluate_models(models,
                     select_model = numTopics,
                     return_model = True,
                     metrics = ['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics = False)
plt.tight_layout()
plt.savefig('model_selection.pdf')  

#Add model to cistopic object.
cistopic_obj_our_groups_only.add_LDA_model(model)

#############
#Visualization
from pycisTopic.clust_vis import run_umap
color_dict_line = {
    'A': '#9A031E',
    'B': '#C75146',
    'C': '#FFA987',
    'D': '#222E50',
    'E': '#8BB174',
}
#Estos son los de python "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"

run_umap(cistopic_obj_our_groups_only, target = 'cell', scale = True)

cistopic_obj_our_groups_only_for_umap = cistopic_obj_our_groups_only

from pycisTopic.clust_vis import plot_metadata
plot_metadata(
    cistopic_obj_our_groups_only,
    reduction_name = 'UMAP',
    color_dictionary = {'group_A': color_dict_line},
    variables = ['group_A'],
    figsize = (10, 10))

plt.savefig('UMAP_cistopic_all_cells.pdf')

#Colorear por d�as (QC)
color_dict_state = {
    'basal': '#9A031E',
    'control_day7': '#f79489',
    'control_day13': '#41729f',
    'control_day20': '#32cd30',
    'ba_day7': '#fac0b9',
    'ba_day13': '#c3e0e5',
    'ba_day20': '#b2d2a4',
}
from pycisTopic.clust_vis import plot_metadata
plot_metadata(
    cistopic_obj_our_groups_only,
    reduction_name = 'UMAP',
    color_dictionary = {'batch': color_dict_state},
    variables = ['batch'],
    figsize = (10, 10))

plt.savefig(work_dir + 'UMAP_cistopic_by_day_QC_all_cells.pdf')


#CARGAMOS LOS PRIMEROS QUE GENERAMOS
#cistopic_obj = pickle.load(open(os.path.join(work_dir + '/Session_objects/cistopic_obj.pkl'), 'rb'))
#region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir + '/Session_objects/region_bin_topics_otsu.pkl'), 'rb'))
#markers_dict = pickle.load(open(os.path.join(work_dir + '/Session_objects/markers_dict.pkl'), 'rb'))

#Inferring candidate enhancer regions
#Next we will infer candidate enhancer regions by binarization of region-topic probabilities ad calculating differentially accessible regions.
#These regions will be used for the next step, pycistarget, in which we will look wich motifs are enriched in these regions.
from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj_our_groups_only, method='otsu')

#Calculate differential accessible regions (DARs).
#We will calculate DARs for each line (i.e. each line vs all other lines), 
#for each state (i.e. each state vs all other states) and for the specific contrast.
#Imputamos/normalizamos/buscamos regiones variables
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj_our_groups_only, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)

#Compute contrasts
#print('Calculating DARs for each group...')
#markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='group_A', var_features=variable_regions, split_pattern = '-', contrasts = contrasts)

print('Calculating DARs for the contrast for C and T for the A-B-C-D-E developmental trajectory')
contrasts = [[['C_group_B'], ['C_group_A']], [['C_group_C'], ['C_group_B']], [['C_group_D'], ['C_group_C']], [['C_group_E'], ['C_group_D']],
             [['T_group_B'], ['T_group_A']], [['T_group_C'], ['T_group_B']], [['T_group_D'], ['T_group_C']], [['T_group_E'], ['T_group_D']]]

markers_dict = find_diff_features(cistopic_obj_our_groups_only, imputed_acc_obj, variable='contrasts_groups', var_features=variable_regions, split_pattern = '-', contrasts = contrasts)

#############
#Motif enrichment analysis using pycistarget

##########################################################################################
#Save session and restore
import pickle    
        
if not os.path.exists(os.path.join(work_dir + '/Session_objects')):
    os.makedirs(os.path.join(work_dir + '/Session_objects'))   

#Export cistopic object             
pickle.dump(cistopic_obj_our_groups_only, open(os.path.join(work_dir + '/Session_objects/cistopic_obj_our_groups_only.pkl'), 'wb'))

#Export cistopic model   
pickle.dump(models,open(os.path.join(work_dir, 'Session_objects/models.pkl'), 'wb'))

#Export otsu binarized
pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir + '/Session_objects/region_bin_topics_otsu.pkl'), 'wb'))
            
#Export cistopic markers
pickle.dump(markers_dict, open(os.path.join(work_dir + '/Session_objects/markers_dict.pkl'), 'wb'))

##
#Load session
cistopic_obj_our_groups_only = pickle.load(open(os.path.join(work_dir + '/Session_objects/cistopic_obj.pkl'), 'rb'))
models = pickle.load(open(os.path.join(work_dir + '/Session_objects/models.pkl'), 'rb'))
region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir + '/Session_objects/region_bin_topics_otsu.pkl'), 'rb'))
markers_dict = pickle.load(open(os.path.join(work_dir + '/Session_objects/markers_dict.pkl'), 'rb'))
##########################################################################################

##############
#Convert to dictionary of pyranges objects.
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['DARs_contrasts'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs_contrasts'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')

#Que base de datos de picos de referencia: 
#For this analysis we will make use a custom made cistarget database on the consensus peaks.
#A database of motif-to-tf annotation database

##############
#download resources
## Download motif database:
#wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather -P /datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes/DB_resources/

#Download scores database
#wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather -P /datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes/DB_resources/

##Download motif-to-tf annotation database
#wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl -P /datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes/DB_resources/

##############
#Load resources
db_path = '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/DB_resources/'
rankings_db = os.path.join(db_path + 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(db_path + 'hg38_screen_v10_clust.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(db_path + 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')

from scenicplus.wrappers.run_pycistarget import run_pycistarget
sys.stderr = open(os.devnull, "w")  # silence stderr
run_pycistarget(
    region_sets = region_sets,
    species = 'homo_sapiens',
    save_path = os.path.join(work_dir + 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 10,
    _temp_dir = os.path.join(tmp_dir + 'ray_spill'),
    annotation_version = 'v10nr_clust')
sys.stderr = sys.__stderr__  # unsilence stderr


##############
##############
#Inferring enhancer-driven Gene Regulatory Networks (eGRNs) using SCENIC+
#We now have completed all the steps for running the SCENIC+ analysis.
#We will start by creating a scenicplus object containing all the analysis we have done up to this point.
import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')
work_dir = '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/'
tmp_dir = '/datos_2/temp/'

##############
#R code
##############
#Generate h5.ad from seurat scRNAseq data
library('SeuratDisk')
setwd("/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/")
load("temporal_2.RData")

#First select ONLY NEURONS cells then  generate 'group_A' on scRNA data metadata
Idents(rna.combined.seurat) <- "seurat_clusters"
neuro_clusters = c("0","1","2","3","4","6","7","8","13","15")
rna.combined.seurat = subset(x = rna.combined.seurat, idents = neuro_clusters)

rna.combined.seurat@meta.data$group_A = rna.combined.seurat@meta.data$seurat_clusters
rna.combined.seurat@meta.data$group_A[rna.combined.seurat@meta.data$group_A == "4"] <- "A"
rna.combined.seurat@meta.data$group_A[rna.combined.seurat@meta.data$group_A == "0"] <- "B"
rna.combined.seurat@meta.data$group_A[rna.combined.seurat@meta.data$group_A == "7"] <- "C"
rna.combined.seurat@meta.data$group_A[rna.combined.seurat@meta.data$group_A == "3"] <- "D"
rna.combined.seurat@meta.data$group_A[rna.combined.seurat@meta.data$group_A == "13"] <- "E"
rna.combined.seurat@meta.data$group_A[rna.combined.seurat@meta.data$group_A == "15"] <- "E"
rna.combined.seurat@meta.data$group_A[rna.combined.seurat@meta.data$group_A == "1"] <- "E"
rna.combined.seurat@meta.data$group_A[rna.combined.seurat@meta.data$group_A == "8"] <- "E"
rna.combined.seurat@meta.data$group_A[rna.combined.seurat@meta.data$group_A == "2"] <- "E"
rna.combined.seurat@meta.data$group_A[rna.combined.seurat@meta.data$group_A == "6"] <- "E"

SaveH5Seurat(rna.combined.seurat, filename = "scRNA_Seurat_to_python.h5Seurat", assays = "RNA")
Convert("scRNA_Seurat_to_python.h5Seurat", dest = "h5ad")
##############
##############

adata = sc.read_h5ad(os.path.join('/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/scRNA_Seurat_to_python.h5ad'))
cistopic_obj = cistopic_obj_our_groups_only  #cistopic_obj = dill.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
menr = dill.load(open(os.path.join(work_dir + 'motifs/menr.pkl'), 'rb'))

#Before we are able to run SCENIC+ we have to combine the scATAC-seq and the scRNA-seq data into a pseudo multiome dataset.
#The way we do this is by randomly sampling a number of cells from the scATAC-seq and scRNA-seq data for each cell type annotation (in this case this is each cell line). 
#We then average the scRNA-seq and scATAc-seq counts within these cells (of the same cell line) to generate metacells containing both scATAC-seq and scRNA-seq data.
from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np
scplus_obj = create_SCENICPLUS_object(
        GEX_anndata = adata,
        cisTopic_obj = cistopic_obj,
        menr = menr,
        multi_ome_mode = False,
        key_to_group_by = 'group_A',
        nr_cells_per_metacells = 5)

#Print all cell_lines in this case STATES are mathching between both omics and will be used
print(f"The cell lines for which we have scRNA-seq data are:\t{', '.join(set(adata.obs['group_A']) - set(['-']))}")
print(f"The cell lines for which we have scATAC-seq data are:\t{', '.join(set(cistopic_obj.cell_data['group_A']))}")
print(f"The cell lines for which we have both:\t{', '.join(set(cistopic_obj.cell_data['group_A']) & set(adata.obs['group_A']))}")

#Now we can run SCENIC+ as usual.
#First let�s check with which biomart host our gene names match.

ensembl_version_dict = {'110': 'http://www.ensembl.org',
                        '109': 'http://feb2023.archive.ensembl.org/',
                        '108': 'http://oct2022.archive.ensembl.org/',
                        '107': 'http://jul2022.archive.ensembl.org/',
                        '106': 'http://apr2022.archive.ensembl.org/',
                        '105': 'http://dec2021.archive.ensembl.org/',
                        '104': 'http://may2021.archive.ensembl.org/',
                        '103': 'http://feb2021.archive.ensembl.org/',
                        '102': 'http://nov2020.archive.ensembl.org/',
                        '101': 'http://aug2020.archive.ensembl.org/',
                        '100': 'http://apr2020.archive.ensembl.org/',
                        '99': 'http://jan2020.archive.ensembl.org/',
                        '98': 'http://sep2019.archive.ensembl.org/',
                        '97': 'http://jul2019.archive.ensembl.org/',
                        '96': 'http://apr2019.archive.ensembl.org/',
                        '95': 'http://jan2019.archive.ensembl.org/',
                        '94': 'http://oct2018.archive.ensembl.org/',
                        '93': 'http://jul2018.archive.ensembl.org/',
                        '92': 'http://apr2018.archive.ensembl.org/',
                        '91': 'http://dec2017.archive.ensembl.org/',
                        '90': 'http://aug2017.archive.ensembl.org/',
                        '89': 'http://may2017.archive.ensembl.org/',
                        '88': 'http://mar2017.archive.ensembl.org/',
                        '87': 'http://dec2016.archive.ensembl.org/',
                        '86': 'http://oct2016.archive.ensembl.org/',
                        '80': 'http://may2015.archive.ensembl.org/',
                        '77': 'http://oct2014.archive.ensembl.org/',
                        '75': 'http://feb2014.archive.ensembl.org/',
                        '54': 'http://may2009.archive.ensembl.org/'}

import pybiomart as pbm
def test_ensembl_host(scplus_obj, host, species):
    dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Chromosome'] = annot['Chromosome'].astype('str')
    filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    gene_names_release = set(annot['Gene'].tolist())
    ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
    print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
    return ov

n_overlap = {}
for version in ensembl_version_dict.keys():
    print(f'host: {version}')
    try:
        n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'hsapiens')
    except:
        print('Host not reachable')

v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")

#Lests select the more overlap
biomart_host = "http://jul2018.archive.ensembl.org"

#Before running we will also download a list of known human TFs from the human transcription factors database.
#wget -O /datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes/DB_resources/utoronto_human_tfs_v_1.01.txt  http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt

#We will also download a the program bedToBigBed this will be used to generate files which can be uploaded to the UCSC genome browser
#wget -O /datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes/bedToBigBed/bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
#chmod +x /datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes/bedToBigBed/bedToBigBed

#NOT DONE YET !! only keep the first two columns of the PCA embedding in order to be able to visualize this in SCope
#scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]
#scplus_obj.dr_cell['GEX_rep'] = scplus_obj.dr_cell['GEX_rep'].iloc[:, 0:2]

#Now we are ready to run the analysis and let�s run the SCENIC+ workflow
from scenicplus.wrappers.run_scenicplus import run_scenicplus

from scenicplus.wrappers.run_scenicplus import run_scenicplus
try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['group_A'],
        species = 'hsapiens',
        assembly = 'hg38',
        tf_file = '/datos_2/Analysis/ALL_DATA_ANALYSIS/final_final/removing_doublets/new_phase_plots/GSEA/PAPER_PLOTS/ATAC/lvl_1_neuro/scenic_outcomes_control/DB_resources/utoronto_human_tfs_v_1.01.txt',
        save_path = os.path.join(work_dir, 'scenicplus_objects'),
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = True,
        export_to_loom_file = True,
        export_to_UCSC_file = True,
        path_bedToBigBed = work_dir,
        n_cpu = 12,
        _temp_dir = os.path.join(tmp_dir, 'ray_spill'))
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)
    raise(e)

############
#Note on the output of SCENIC+
############

#Both the raw gene expression counts and chromatin accessibility data are stored in
scplus_obj.to_df('EXP').head()
scplus_obj.to_df('ACC').head()

#Cell metatdata is stored in
scplus_obj.metadata_cell.head()

#Region metadata is stored in
scplus_obj.metadata_regions.head()

#Gene metadata is stored in
scplus_obj.metadata_genes.head()

#Motif enrichment data is stored in
scplus_obj.menr.keys()

#Dimensionality reductions of the cells are stored in
scplus_obj.dr_cell.keys()

#Additional unstructured data will be stored in

#Pre-SCENIC+
#Cistromes: this contains TFs together with target regions based on the motif enrichment analysis (i.e. prior to running SCENIC+)
#search_space: this is a dataframe containing the search space for each gene.
#region_to_gene: this is a dataframe containing region to gene links prior to running SCENIC+ (i.e unfiltered/raw region to gene importance scores and correlation coefficients).
#TF2G_adj: this is a datafram containing TF to gene links prior to running SCENIC+ (i.e unfiltered/raw TF to gene importance scores and correlation coefficients).

#POST-SCENIC+
#eRegulons: this is the raw output from the SCENIC+ analysis. We will go into a bit more detail for these below.
#eRegulon_metadata: this is a dataframe containing the same information as eRegulons bit in a format which is a bit easier to parse for a human.
#eRegulon_signatures: this is a dictionary with target regions and genes for each eRegulon
#eRegulon_AUC: this slot contains dataframes with eRegulon enrichment scores calculated using AUCell (see below).
#pseudobulk: contains pseudobulked gene expression and chromatin accessibility data, this is used to calculated TF to eRegulon correlation values.
#TF_cistrome_correlation: contains correlation values between TF expression and eRegulon enrichment scores (seperate entries for target gene and target region based scores).
#eRegulon_AUC_thresholds: contains thresholds on the AUC values (eRegulon enrichment scores), this is necessary to be able to visualize the results in SCope
#RSS: contains eRegulon Specificity Scores (RSS), a measure on how cell type specific an eRegulon is.
#DEGs: contains Differentially Expressed Genes.
#DARs: contains Differentially Accessibile Regions.
scplus_obj.uns.keys()

#The main output of SCENIC+ are eRegulons
#This is initially stored in a list of eRegulon classes as depicted below
scplus_obj.uns['eRegulons'][0:5]

#each eRegulon has the following information (attributes):
#cistrome_name: name of the cistrome (from scenicplus.uns['Cistromes']) from which this eRegulon was created.
#context: specifies the binarization method(s) used for binarizing region to gene relationships and wether positive/negative region-to-gene and TF-to-gene relationships were used
#gsea_adj_pval/gsea_enrichment_score/gsea_pval/in_leading_edge: are internal parameters used for generating the eRegulons. The values are lost when generating the final eRegulons because results from several analysis (different binarization methods) are combined.
#is_extended: specifies wether extended (i.e. non-direct) motif-to-TF annotations were used.
#n_target_genes: number of target genes.
#n_target_regions: number of target regions.
#regions2genes: region to gene links after running SCENIC+.
#target_genes: target genes of the eRegulon
#target_regions: target regions of the eRegulon
#transcription_factor: TF name
for attr in dir(scplus_obj.uns['eRegulons'][0]):
    if not attr.startswith('_'):
        print(f"{attr}: {getattr(scplus_obj.uns['eRegulons'][0], attr) if not type(getattr(scplus_obj.uns['eRegulons'][0], attr)) == list else getattr(scplus_obj.uns['eRegulons'][0], attr)[0:5]}")

#The information of all eRegulons is combined in the eRegulon_metadata dataframe
scplus_obj.uns['eRegulon_metadata'].head()

#For the eRegulon names we use the following convetion:
#<TF NAME>_<TF-TO-GENE RELATIONSHIP (+/-)>_<REGION-TO-GENE RELATIONSHIP (+/-)>_(NUMBER OF TARGET REGIONS(r)/GENES(g))
#For example the name: ARID3A_+_+_(364r) and ARID3A_+_+_(278g) indicates that the we found an eRegulon for the TF ARID3A which has 364 target regions and 278 target genes, that expression of the TF correlates positively with the expression of all the target genes (first + sign) and that the accessibility of all target regions correlates positively with the expression of all target regions (seconf + sign).

###########
#Downstream analysis
###########
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

import dill
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus_objects/scplus_obj.pkl'), 'rb'))

#Simplifying and filtering SCENIC+ output
#Given the multitude of eRegulons that can be generated for each TF (see above) we will first simplify the result by: 
#1. Only keeping eRegulons with an extended annotation if there is no direct annotation available (given that the confidence of direct motif annotations is in genral higher). 
#2. Discarding eRegulons for which the region-to-gene correlation is negative (these are often noisy). 
#3. Renaming the eRegulons so that eRegulons with the suffix TF_+_+ become TF_+ and those with TF_-_+ become TF_-.
from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
apply_std_filtering_to_eRegulons(scplus_obj)

#This will create two new entries in the scenicplus object: 
#scplus_obj.uns['eRegulon_metadata_filtered']  
#scplus_obj.uns['eRegulon_signatures_filtered'] containing the simplified results. We will use these for downstream analysis.
scplus_obj.uns['eRegulon_metadata_filtered'].head()

#eRegulon enrichment scores
#We can score the enrichment of eRegulons using the AUCell function
#This function takes as input a gene or region based ranking (ranking of genes/regions based on the expression/accessibility per cell) and a list of eRegulons.
#These values were already calculated in the wrapper function but let�s recalculate them using the filtered output.
from scenicplus.eregulon_enrichment import score_eRegulons
region_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus_objects/region_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function
gene_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus_objects/gene_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function
score_eRegulons(scplus_obj,
                ranking = region_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type= 'region',
                auc_threshold = 0.05,
                normalize = False,
                n_cpu = 5)
score_eRegulons(scplus_obj,
                gene_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type = 'gene',
                auc_threshold = 0.05,
                normalize= False,
                n_cpu = 5)

#eRegulon dimensionality reduction
#Based on the enrichment scores calculated above we can generate dimensionality reductions (e.g. tSNE and UMAP).
#To calculate these dimensionality reductions we use both the regions and gene based enrichment scores.
from scenicplus.dimensionality_reduction import run_eRegulons_tsne, run_eRegulons_umap
run_eRegulons_umap(
    scplus_obj = scplus_obj,
    auc_key = 'eRegulon_AUC_filtered',
    reduction_name = 'eRegulons_UMAP', #overwrite previously calculated UMAP
)
run_eRegulons_tsne(
    scplus_obj = scplus_obj,
    auc_key = 'eRegulon_AUC_filtered',
    reduction_name = 'eRegulons_tSNE', #overwrite previously calculated tSNE
)

#Let�s visualize the UMAP and tSNE stored respectively in eRegulons_UMAP and eRegulons_tSNE, these are calculated based on the combined region and gene AUC values described above.
#Let�s also add some nice colours by specifying a color_dictionary.
from scenicplus.dimensionality_reduction import plot_metadata_given_ax
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline

#specify color_dictionary
color_dict = {
    'A': "#065143",
    'B': "#70B77E",
    'C': "#E0A890",
    'D': "#F56476",
    'E': "#CE1483"
}
color_dict_line = {
    'A': '#9A031E',
    'B': '#C75146',
    'C': '#FFA987',
    'D': '#222E50',
    'E': '#8BB174',
}

full_name = work_dir + 'SCENIC+' + '_eRegulon_UMAP_and_TSNE' + '.pdf'
fig, axs = plt.subplots(ncols=2, figsize = (16, 8))
plot_metadata_given_ax(
    scplus_obj=scplus_obj,
    ax = axs[0],
    reduction_name = 'eRegulons_UMAP',
    variable = 'group_A', #note the GEX_ prefix, this metadata originated from the gene expression metadata (on which we did the cell type annotation before)
    color_dictionary={'GEX_celltype': color_dict}
)
plot_metadata_given_ax(
    scplus_obj=scplus_obj,
    ax = axs[1],
    reduction_name = 'eRegulons_tSNE',
    variable = 'group_A', #note the GEX_ prefix, this metadata originated from the gene expression metadata (on which we did the cell type annotation before)
    color_dictionary={'GEX_celltype': color_dict}
)
fig.tight_layout()
sns.despine(ax = axs[0]) #remove top and right edge of axis border
sns.despine(ax = axs[1]) #remove top and right edge of axis border
#plt.show()
plt.savefig(full_name)
plt.close()

#plot the activity / expression of an eRegulon on the dimensionality reduction
#Nex we visualize the gene expression and target gene and region activity of some eRegulons on the tSNE.
scplus_obj.uns['eRegulon_metadata'].head()

from scenicplus.dimensionality_reduction import plot_eRegulon

full_name = work_dir + 'SCENIC+' + '_eRegulon_gene_geneExp_region_activity__Inpairment_1' + '.pdf'
plot_eRegulon(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_tSNE',
    selected_regulons = ['EBF1_+','TCF7L1_+','FOS_+','ASCL1_+'],
    scale = True,
    auc_key = 'eRegulon_AUC_filtered')
plt.savefig(full_name)
plt.close()

full_name = work_dir + 'SCENIC+' + '_eRegulon_gene_geneExp_region_activity__Inpairment_2' + '.pdf'
plot_eRegulon(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_tSNE',
    selected_regulons = ['POU2F2_+','PAX6_+','SOX2_+','PBX3_+','EGR1_+','GATA2_+','ONECUT1_+','ONECUT2_+','MEIS1_+','NKX6-1_+'],
    scale = True,
    auc_key = 'eRegulon_AUC_filtered')
plt.savefig(full_name)
plt.close()

full_name = work_dir + 'SCENIC+' + '_eRegulon_gene_geneExp_region_activity__Inpairment_3' + '.pdf'
plot_eRegulon(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_tSNE',
    selected_regulons = ['GATA2_+','ONECUT1_+','ONECUT2_+','MEIS1_+','NKX6-1_+'],
    scale = True,
    auc_key = 'eRegulon_AUC_filtered')
plt.savefig(full_name)
plt.close()

full_name = work_dir + 'SCENIC+' + '_eRegulon_gene_geneExp_region_activity__Inpairment_Mature' + '.pdf'
plot_eRegulon(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_tSNE',
    selected_regulons = ['PAX6_+'],
    scale = True,
    auc_key = 'eRegulon_AUC_filtered')
plt.savefig(full_name)
plt.close()

#We can also plot only the activity of an eRegulon
fig, ax = plt.subplots(figsize = (8,8))
plot_AUC_given_ax(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_tSNE',
    feature = 'PAX5_+_(119g)',
    ax = ax,
    auc_key = 'eRegulon_AUC_filtered',
    signature_key = 'Gene_based')
sns.despine(ax = ax)
plt.show()


#For eRegulons it is often usefull to visualize both information on the TF/target genes expression and region accessibility at the same time.
#dotplot-heatmap
#A dotplot-heatmap is a useful way to visualize this. Here the color of the heatmap can be used to visualize one aspect of the eRegulon (for example TF expression) and the size of the dot can be used to visualize another aspect (for example the enrichment (AUC value) of eRegulon target regions).
#Before we plot the the dotplot-heatmap let�s first select some high quality eRegulons to limit the amount of space we need for the plot. One metric which can be used for selecting eRegulons is the correlation between TF expression and target region enrichment scores (AUC values). Let�s (re)calculate this value based on the simplified eRegulons
#We first generate pseudobulk gene expression and region accessibility data, per celltype, to limit the amount of noise for the correlation calculation.
from scenicplus.cistromes import TF_cistrome_correlation, generate_pseudobulks

generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'group_A',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Gene_based')
        
generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'group_A',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Region_based')

TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'group_A',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Gene_based',
            out_key = 'filtered_gene_based')
            
TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'group_A',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Region_based',
            out_key = 'filtered_region_based')

scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].head()

#Let�s visualize these correlations in a scatter plot and select eRegulons for which the correlaiton coefficient is above 0.70 or below -0.75
import numpy as np
n_targets = [int(x.split('(')[1].replace('r)', '')) for x in scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Cistrome']]
rho = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'].to_list()
adj_pval = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Adjusted_p-value'].to_list()

thresholds = {
        'rho': [-0.75, 0.70],
        'n_targets': 0
}
import seaborn as sns

full_name = work_dir + 'SCENIC+' + '_eRegulon_Correlation_+0.70_antiCorrelation_-0.75' + '.pdf'
fig, ax = plt.subplots(figsize = (10, 5))
sc = ax.scatter(rho, n_targets, c = -np.log10(adj_pval), s = 5)
ax.set_xlabel('Correlation coefficient')
ax.set_ylabel('nr. target regions')
#ax.hlines(y = thresholds['n_targets'], xmin = min(rho), xmax = max(rho), color = 'black', ls = 'dashed', lw = 1)
ax.vlines(x = thresholds['rho'], ymin = 0, ymax = max(n_targets), color = 'black', ls = 'dashed', lw = 1)
ax.text(x = thresholds['rho'][0], y = max(n_targets), s = str(thresholds['rho'][0]))
ax.text(x = thresholds['rho'][1], y = max(n_targets), s = str(thresholds['rho'][1]))
sns.despine(ax = ax)
fig.colorbar(sc, label = '-log10(adjusted_pvalue)', ax = ax)
plt.savefig(full_name)
plt.close()

#Select eRegulons base on correlation
selected_cistromes = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].loc[
        np.logical_or(
                scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] > thresholds['rho'][1],
                scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] < thresholds['rho'][0]
        )]['Cistrome'].to_list()
selected_eRegulons = [x.split('_(')[0] for x in selected_cistromes]
selected_eRegulons_gene_sig = [
        x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()
        if x.split('_(')[0] in selected_eRegulons]
selected_eRegulons_region_sig = [
        x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()
        if x.split('_(')[0] in selected_eRegulons]
#save the results in the scenicplus object
scplus_obj.uns['selected_eRegulon'] = {'Gene_based': selected_eRegulons_gene_sig, 'Region_based': selected_eRegulons_region_sig}
print(f'selected: {len(selected_eRegulons_gene_sig)} eRegulons')
###########
#66 eRegulons selected base ion this filters
###########

#Let�s save these changes we have made to the scenicplus_obj
dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus_objects/scplus_obj_after_filtering.pkl'), 'wb'), protocol=-1)

#Let�s plot the heatmap-dotplot for the selected most interesting eRegulons
from scenicplus.plotting.dotplot import heatmap_dotplot

full_name = work_dir + 'SCENIC+' + '_Heatmap_of_filtered_most_relevant_eRegulons' + '.pdf'
heatmap_dotplot(
        scplus_obj = scplus_obj,
        size_matrix = scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'], #specify what to plot as dot sizes, target region enrichment in this case
        color_matrix = scplus_obj.to_df('EXP'), #specify  what to plot as colors, TF expression in this case
        scale_size_matrix = True,
        scale_color_matrix = True,
        group_variable = 'group_A',
        subset_eRegulons = scplus_obj.uns['selected_eRegulon']['Gene_based'],
        index_order = ['A', 'B', 'C', 'D', 'E'],
        figsize = (5, 20),
        orientation = 'vertical')
plt.savefig(full_name)
plt.close()

#Overlap of predicted target regions
#An interesting aspect of gene regulation is transcription factor cooperativity (i.e. multiple TFs cobinding the same enhancer together driving gene expression).
#By looking at the overlap of predicted target regions of TFs we can infer potential cooperativity events.
#Let�s look at the overlap of target regions of the top 5 TFs per cell type based on the Regulon Specificity Score (RSS).
#First we calculate the RSS for the target regions of the selected eRegulons.
from scenicplus.RSS import *
regulon_specificity_scores(
        scplus_obj,
        variable = 'group_A',
        auc_key = 'eRegulon_AUC_filtered',
        signature_keys = ['Region_based'],
        selected_regulons = [x for x in scplus_obj.uns['selected_eRegulon']['Region_based'] if '-' not in x],
        out_key_suffix = '_filtered')

#Let�s visualize the RSS values using a scatter plot
full_name = work_dir + 'SCENIC+' + '_transcription_factor_cooperativity_RSS' + '.pdf'
plot_rss(scplus_obj, 'group_A_filtered', num_columns=5, top_n=10, figsize = (60, 12))
plt.savefig(full_name)
plt.close()

#Next we select the top 10 eRegulons per cell type (group_A)
flat_list = lambda t: [item for sublist in t for item in sublist]
selected_markers = list(set(flat_list(
    [scplus_obj.uns['RSS']['group_A_filtered'].loc[celltype].sort_values(ascending = False).head(10).index.to_list()
    for celltype in scplus_obj.uns['RSS']['group_A_filtered'].index])))
    
from scenicplus.plotting.correlation_plot import *
region_intersetc_data, Z = jaccard_heatmap(
        scplus_obj,
        method = 'intersect',
        gene_or_region_based = 'Region_based',
        use_plotly = True,
        selected_regulons = selected_markers,
        signature_key = 'eRegulon_signatures_filtered',
        figsize = (10, 10), return_data = True, vmax = 0.5, cmap = 'plasma')    
  
#Plot top 10 eRegulons per cell type 
import seaborn as sns
full_name = work_dir + 'SCENIC+' + '_top10_RSS_eRegulons_per_group' + '.pdf' 
sns.heatmap(region_intersetc_data, cmap='plasma')
plt.savefig(full_name)
plt.close()  

#Plotting a network
#eRegulons can also be visualized in a network. Simple plots can be made using python. For more complicated plots (i.e. containing many nodes and edges) we suggest exporting your network to cytoscape.
#Let�s create a very simple network for a cell type. We will use the top 1000 highly variable regions and genes in this plot. If you want to use more feautures please export your nework to cytoscape.
from pycisTopic.diff_features import find_highly_variable_features
hvr = find_highly_variable_features(scplus_obj.to_df('ACC').loc[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Region']))], n_top_features=1000, plot = False)
hvg = find_highly_variable_features(scplus_obj.to_df('EXP')[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Gene']))].T, n_top_features=1000, plot = False)

#First we format the eRegulons into a table which can be used to create a network using the package 'networkx'
from scenicplus.networks import create_nx_tables, create_nx_graph, plot_networkx, export_to_cytoscape
nx_tables = create_nx_tables(
    scplus_obj = scplus_obj,
    eRegulon_metadata_key ='eRegulon_metadata_filtered',
    subset_eRegulons = ['PAX6', 'EBF1', 'SOX2'],
    subset_regions = hvr,
    subset_genes = hvg,
    add_differential_gene_expression = True,
    add_differential_region_accessibility = True,
    differential_variable = ['group_A'])

#Next we layout the graph.
G, pos, edge_tables, node_tables = create_nx_graph(nx_tables,
                   use_edge_tables = ['TF2R','R2G'],
                   color_edge_by = {'TF2R': {'variable' : 'TF', 'category_color' : {'PAX6': 'Orange', 'EBF1': 'Purple', 'SOX2': 'Red'}},
                                    'R2G': {'variable' : 'R2G_rho', 'continuous_color' : 'viridis', 'v_min': -1, 'v_max': 1}},
                   transparency_edge_by =  {'R2G': {'variable' : 'R2G_importance', 'min_alpha': 0.1, 'v_min': 0}},
                   width_edge_by = {'R2G': {'variable' : 'R2G_importance', 'max_size' :  1.5, 'min_size' : 1}},
                   color_node_by = {'TF': {'variable': 'TF', 'category_color' : {'PAX6': 'Orange', 'EBF1': 'Purple', 'SOX2': 'Red'}},
                                    'Gene': {'variable': 'group_A_Log2FC_A', 'continuous_color' : 'bwr'},
                                    'Region': {'variable': 'group_A_Log2FC_A', 'continuous_color' : 'viridis'}},
                   transparency_node_by =  {'Region': {'variable' : 'group_A_Log2FC_A', 'min_alpha': 0.1},
                                            'Gene': {'variable' : 'group_A_Log2FC_A', 'min_alpha': 0.1}},
                   size_node_by = {'TF': {'variable': 'fixed_size', 'fixed_size': 30},
                                    'Gene': {'variable': 'fixed_size', 'fixed_size': 15},
                                    'Region': {'variable': 'fixed_size', 'fixed_size': 10}},
                   shape_node_by = {'TF': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
                                    'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
                                    'Region': {'variable': 'fixed_shape', 'fixed_shape': 'diamond'}},
                   label_size_by = {'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 15.0},
                                    'Gene': {'variable': 'fixed_label_size', 'fixed_label_size': 5.0},
                                    'Region': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0}},
                   layout='kamada_kawai_layout',
                   scale_position_by=250)

#Finally we can visualize the network.
#In this network diamond shapes represent regions and they are color coded by their log2fc value in a cell type target genes and TFs are visualized using circles and are labeled.
full_name = work_dir + 'SCENIC+' + '_eRegulons_network_group__A' + '.pdf'
plt.figure(figsize=(10,10))
plot_networkx(G, pos)
plt.savefig(full_name)
plt.close()  

#We can also export this network to a format which can be opened in Cytoscape.
export_to_cytoscape(G, pos, out_file = os.path.join(work_dir, 'scenicplus_objects/network_B_cells.cys'))
















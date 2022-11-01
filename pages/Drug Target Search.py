import streamlit as st
import pandas as pd

def convert_df(df):
     # IMPORTANT: Cache the conversion to prevent computation on every rerun
     return df.to_csv(index = False)

def drug_search(genes_list, df):
    
    # check gene list to see if user genes are in druggable genome
    no_gene = []
    yes_gene = []
    for gene in genes_list:
        if gene not in list(df['hgnc_names']):
            no_gene.append(gene)
        else:
            yes_gene.append(gene)
            
    # let user know if any of their genes are not drug targets     
    if len(no_gene) == 1:
        no_result = f'{no_gene[0]} is not a known therapeutic target. Please check your target gene for spelling/formatting and try again.'
    elif len(no_gene) > 1:
        no_result = f'{", ".join(no_gene)} are not known therapeutics targets. Please check your target genes for spelling/formatting and try again.'
    else:
        no_result = None
    
    # return known drug targets
    results = df.query("hgnc_names == @yes_gene")
    print(f'{", ".join(yes_gene)} are known drug targets!')
    
    return results, no_result


# session state variables for  gene dataframe
if 'genes_list' not in st.session_state: # user provided genes
    st.session_state['genes_list'] = None
if 'gene_results_df' not in st.session_state:
    st.session_state['gene_results_df'] = None
if 'gene_results_run' not in st.session_state:
    st.session_state['gene_results_run'] = None

# session state variables for druggable gene dataframe
if 'drugdf' not in st.session_state:
    st.session_state['drugdf'] = None


# load in druggable genome list
drug_df = pd.read_csv('./druggable_genome.csv', sep = ',') 
st.session_state['drugdf'] = drug_df

st.title('Therapeutic Gene Target Search')
st.write('Provide a gene or list of genes to search against our therapeutic target database.')

# user input - genes list
genes_list = st.text_area('Input a gene or list of genes (seperated by comma):')
# split list
genes_list = genes_list.split(',')
new_genes = []
for gene in genes_list:
    new_genes.append(gene.strip())

st.session_state['genes_list'] = new_genes # add to session state

# user submits list
search_run = st.button('Submit')

if search_run:
    st.session_state['gene_results_run'] = 'run'

if st.session_state['gene_results_run'] == 'run':
    gene_results, no_gene = drug_search(st.session_state['genes_list'], st.session_state['drugdf'])

    if no_gene != None:
        st.write(no_gene)
        st.session_state['gene_results_df'] = gene_results
    else:
        st.session_state['gene_results_df'] = gene_results

    # display results dataframe
    st.dataframe(st.session_state['gene_results_df'])

    # allow user to download results
    st.download_button(label="Download results as CSV", data=convert_df(st.session_state['gene_results_df']), mime='text/csv')
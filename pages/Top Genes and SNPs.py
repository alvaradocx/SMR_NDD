import streamlit as st
import pandas as pd

def convert_df(df):
     # IMPORTANT: Cache the conversion to prevent computation on every rerun
     return df.to_csv(index = False)

def create_df(df, diseases, omics):
  # determine max possible options 
  dx_max = len(df['Disease'])
  o_max = len(df['Omic'].unique())

  # num of options in provided lists
  dx_len = len(diseases)
  o_len = len(omics)

  if dx_len == dx_max and o_max == o_len:
    return df

  elif dx_len == dx_max:
    new_df = df[df['Omic'].isin(omics)]
    return new_df

  elif o_max == o_len:
    new_df = df[df['Disease'].isin(diseases)]
    return new_df
  
  else:
    new_df = df.query("Disease in @diseases and Omic in @omics")
    return new_df

def fdr(df):
    
    # pull p_vals
    p_vals = df['p_SMR_multi']
    
    # run fdr on pvals
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    
    # add column to df
    df['FDR_pval'] = fdr

    return df

def get_top_genes(df, diseases, omics_list = None, HLA = False, create_chart = False, heidi = 0.01, pval = 0.05, peqtl = 0.05):
    # preset list of NDD omics
    ndd_omics = ['Brain_Amygdala', 'Brain_Hippocampus', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Hypothalamus',
    'Brain_Cerebellar_Hemisphere', 'Brain_Substantia_nigra', 'Brain_Caudate_basal_ganglia', 'Brain_Putamen_basal_ganglia', 'Brain_Cerebellum',  'Brain_Spinal_cord_cervical_c-1',
    'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Liver', 'Nerve_Tibial', 'Whole_Blood', 'Muscle_Skeletal', 'Cerebellum_metaBrain', 'Spinalcord_metaBrain', 'brain_eMeta', 'Cortex_metaBrain',
    'Basalganglia_metaBrain',  'Hippocampus_metaBrain', 'blood_eQTLgen', 'brain_mMeta', 'blood_Bryois', 'blood_mcrae']
    all_omics = list(df['Omic'].unique())

    if HLA:
         df = df[~df['annotated_gene'].str.contains("HLA")]

    smr_df = pd.DataFrame() 

    # filter out for selected diseases
    for dx in diseases:
        dx_df = df.query(f"Disease == '{dx}'")
        smr_df = pd.concat([smr_df, dx_df])

    # final df
    top_genes_df = pd.DataFrame()

    # apply omics filter from user chosen omics
    if omics_list == 'all':
        for disease in diseases:
            for omic in all_omics:
                temp_df = smr_df.query(f"Disease == '{disease}' & Omic == '{omic}' & p_SMR_multi < {pval} & p_HEIDI > {heidi} & p_eQTL < {peqtl}")
                temp_df = fdr(temp_df)
                top_gene = temp_df[temp_df['FDR_pval'] == temp_df['FDR_pval'].min()]
                top_gene = top_gene[['Omic', 'Disease', 'annotated_gene', 'topRSID', 'b_GWAS','p_GWAS', 'b_SMR', 'p_SMR_multi', 'p_HEIDI', 'p_eQTL', 'FDR_pval']]
                
                top_genes_df = pd.concat([top_genes_df,top_gene])
        
        if create_chart: # return chart
            top_genes_chart = pd.DataFrame(index = all_omics, columns = diseases)
            for index, row in top_genes_df.iterrows():
                omic_ind = row['Omic']
                disease_col = row['Disease']
                if row['annotated_gene'] != 'missing_gene':
                    top_genes_chart.loc[omic_ind, disease_col] = row['annotated_gene']
                else:
                    top_genes_chart.loc[omic_ind, disease_col] = row['topRSID']

            return top_genes_df,top_genes_chart

        else: 
            return top_genes_df
    
    elif omics_list != None: # Custom omics list
        for disease in diseases:
            for omic in omics_list:
                temp_df = smr_df.query(f"Disease == '{disease}' & Omic == '{omic}' & p_SMR_multi < {pval} & p_HEIDI > {heidi} & p_eQTL < {peqtl}")
                temp_df = fdr(temp_df)
                top_gene = temp_df[temp_df['FDR_pval'] == temp_df['FDR_pval'].min()]
                top_gene = top_gene[['Omic', 'Disease', 'annotated_gene', 'topRSID', 'b_GWAS','p_GWAS', 'b_SMR', 'p_SMR_multi', 'p_HEIDI', 'p_eQTL', 'FDR_pval']]
                
                top_genes_df = pd.concat([top_genes_df,top_gene])
        
        if create_chart: # return chart
            top_genes_chart = pd.DataFrame(index = omics_list, columns = diseases)
            for index, row in top_genes_df.iterrows():
                omic_ind = row['Omic']
                disease_col = row['Disease']
                if row['annotated_gene'] != 'missing_gene':
                    top_genes_chart.loc[omic_ind, disease_col] = row['annotated_gene']
                else:
                    top_genes_chart.loc[omic_ind, disease_col] = row['topRSID']

            return top_genes_df,top_genes_chart

        else: 
            return top_genes_df

    else: # if using default ndd_omics listthan no omics_list is provided to the function; acts as default option
        for disease in diseases:
            for omic in ndd_omics:
                temp_df = smr_df.query(f"Disease == '{disease}' & Omic == '{omic}' & p_SMR_multi < {pval} & p_HEIDI > {heidi} & p_eQTL < {peqtl}")
                temp_df = fdr(temp_df)
                top_gene = temp_df[temp_df['FDR_pval'] == temp_df['FDR_pval'].min()]
                top_gene = top_gene[['Omic', 'Disease', 'annotated_gene', 'topRSID', 'b_GWAS','p_GWAS', 'b_SMR', 'p_SMR_multi', 'p_HEIDI', 'p_eQTL', 'FDR_pval']]
                
                top_genes_df = pd.concat([top_genes_df,top_gene])

        if create_chart: # return chart
            top_genes_chart = pd.DataFrame(index = ndd_omics, columns = diseases)
            for index, row in top_genes_df.iterrows():
                omic_ind = row['Omic']
                disease_col = row['Disease']
                if row['annotated_gene'] != 'missing_gene':
                    top_genes_chart.loc[omic_ind, disease_col] = row['annotated_gene']
                else:
                    top_genes_chart.loc[omic_ind, disease_col] = row['topRSID']
            return top_genes_df, top_genes_chart

        else:
            return top_genes_df
    

def get_top_snps(df, diseases, omics_list = None, HLA = False, create_chart = False, heidi = 0.01, pval = 0.05, peqtl = 0.05):
     # preset list of NDD omics
    ndd_omics = ['Brain_Amygdala', 'Brain_Hippocampus', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Hypothalamus',
    'Brain_Cerebellar_Hemisphere', 'Brain_Substantia_nigra', 'Brain_Caudate_basal_ganglia', 'Brain_Putamen_basal_ganglia', 'Brain_Cerebellum',  'Brain_Spinal_cord_cervical_c-1',
    'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Liver', 'Nerve_Tibial', 'Whole_Blood', 'Muscle_Skeletal', 'Cerebellum_metaBrain', 'Spinalcord_metaBrain', 'brain_eMeta', 'Cortex_metaBrain',
    'Basalganglia_metaBrain',  'Hippocampus_metaBrain', 'blood_eQTLgen', 'brain_mMeta', 'blood_Bryois', 'blood_mcrae']
    all_omics = list(df['Omic'].unique())

    # remove all genes with HLA in name if user choses
    if HLA:
         df = df[~df['annotated_gene'].str.contains("HLA")]
    

    smr_df = pd.DataFrame() 

    # filter out for selected diseases
    for dx in diseases:
        dx_df = df.query(f"Disease == '{dx}'")
        smr_df = pd.concat([smr_df, dx_df])

    # final df
    top_snps_df = pd.DataFrame()

    # apply omics filter from user chosen omics
    if omics_list == 'all':
        for disease in diseases:
            for omic in all_omics:
                temp_df = smr_df.query(f"Disease == '{disease}' & Omic == '{omic}' & p_SMR_multi < {pval} & p_HEIDI > {heidi} & p_eQTL < {peqtl}")
                temp_df = fdr(temp_df)
                top_snp = temp_df[temp_df['FDR_pval'] == temp_df['FDR_pval'].min()]
                top_snp = top_snp[['Omic', 'Disease', 'annotated_gene', 'topRSID', 'b_GWAS','p_GWAS', 'b_SMR', 'p_SMR_multi', 'p_HEIDI', 'p_eQTL', 'FDR_pval']]
                
                top_snps_df = pd.concat([top_snps_df,top_snp])
        
        if create_chart: # return chart
            top_snps_chart = pd.DataFrame(index = all_omics, columns = diseases)
            for index, row in top_snps_df.iterrows():
                omic_ind = row['Omic']
                disease_col = row['Disease']
                top_snps_chart.loc[omic_ind, disease_col] = row['topRSID']

            return top_snps_df,top_snps_chart

        else: 
            return top_snps_df
    
    elif omics_list != None: # Custom omics list
        for disease in diseases:
            for omic in omics_list:
                temp_df = smr_df.query(f"Disease == '{disease}' & Omic == '{omic}' & p_SMR_multi < {pval} & p_HEIDI > {heidi} & p_eQTL < {peqtl}")
                temp_df = fdr(temp_df)
                top_snp = temp_df[temp_df['FDR_pval'] == temp_df['FDR_pval'].min()]
                top_snp = top_snp[['Omic', 'Disease', 'annotated_gene', 'topRSID', 'b_GWAS','p_GWAS', 'b_SMR', 'p_SMR_multi', 'p_HEIDI', 'p_eQTL', 'FDR_pval']]
                
                top_snps_df = pd.concat([top_snps_df,top_snp])
        
        if create_chart: # return chart
            top_snps_chart = pd.DataFrame(index = all_omics, columns = diseases)
            for index, row in top_snps_df.iterrows():
                omic_ind = row['Omic']
                disease_col = row['Disease']
                top_snps_chart.loc[omic_ind, disease_col] = row['topRSID']

            return top_snps_df,top_snps_chart

        else: 
            return top_snps_df

    else: # if using default ndd_omics list
        for disease in diseases:
            for omic in ndd_omics:
                temp_df = smr_df.query(f"Disease == '{disease}' & Omic == '{omic}' & p_SMR_multi < {pval} & p_HEIDI > {heidi} & p_eQTL < {peqtl}")
                temp_df = fdr(temp_df)
                top_snp = temp_df[temp_df['FDR_pval'] == temp_df['FDR_pval'].min()]
                top_snp = top_snp[['Omic', 'Disease', 'annotated_gene', 'topRSID', 'b_GWAS','p_GWAS', 'b_SMR', 'p_SMR_multi', 'p_HEIDI', 'p_eQTL', 'FDR_pval']]
                
                top_snps_df = pd.concat([top_snps_df,top_snp])

        if create_chart: # return chart
            top_snps_chart = pd.DataFrame(index = all_omics, columns = diseases)
            for index, row in top_snps_df.iterrows():
                omic_ind = row['Omic']
                disease_col = row['Disease']
                top_snps_chart.loc[omic_ind, disease_col] = row['topRSID']

            return top_snps_df,top_snps_chart

        else: 
            return top_snps_df
    
def gene_count(df): # count how many times a gene shows up in the top genes df
    df_counts = df.stack().value_counts().rename_axis('value').reset_index(name='count')

    return df_counts

def drug_identify(user_df, drug_df):

    # pull out unique genes from user df
    all_top_gene = list(user_df['value'])
    all_top_gene_unique = []
    for gene in all_top_gene:
        if gene not in all_top_gene_unique:
            all_top_gene_unique.append(gene)
        else:
            continue
    

    gene_df = pd.DataFrame(all_top_gene_unique, columns = ['gene'])

    drug_merge = drug_df.merge(gene_df, left_on = 'hgnc_names', right_on = 'gene')
    return drug_merge


if 'top_df' not in st.session_state:
    st.session_state['top_df'] = None
if 'top_chart' not in st.session_state:
    st.session_state['top_chart'] = None
if 'query_submit' not in st.session_state:
    st.session_state['query_submit'] = 'not_run'
if 'top_submit' not in st.session_state:
    st.session_state['top_submit'] = 'not_run'
if 'dx_list' not in st.session_state:
    st.session_state['dx_list'] = None
if 'peqtl' not in st.session_state:
    st.session_state['peqtl'] = 0.05
if 'pval' not in st.session_state:
    st.session_state['pval'] = 0.05
if 'heidi' not in st.session_state:
    st.session_state['heidi'] = 0.01
if 'omic_list' not in st.session_state:
    st.session_state['omic_list'] = ''

if 'df_name' not in st.session_state:
    st.session_state['df_name'] = ''
if 'chart_name' not in st.session_state:
    st.session_state['chart_name'] = ''

# session state variables for druggable gene dataframe
if 'drugdf' not in st.session_state:
    st.session_state['drugdf'] = None

# session state variables for value_count + gene druggability
if 'drug_status' not in st.session_state:
    st.session_state['drug_status'] = 'not_run'
if 'value_df' not in st.session_state:
    st.session_state['value_df'] = None
if 'value_status' not in st.session_state:
    st.session_state['value_status'] = 'not_run'
if 'gene_drugdf' not in st.session_state:
    st.session_state['gene_drugdf'] = None
if 'drug_name' not in st.session_state:
    st.session_state['drug_name'] = ''

main_df = st.session_state['main_data']
# load in druggable genome list
drug_df = pd.read_csv('./druggable_genome.csv', sep = ',') 
st.session_state['drugdf'] = drug_df

st.title('Top Genes and Top SNPs')
st.write('This interactive tool will allow you to extract the top genes and/or SNPs across available diseases and available omics. A top candidate is defined as a having the most significant p-value after multiple test correction.')


with st.form('query selection'):
    gene_or_not = st.radio("Extract top genes or top SNPs?", ('Genes', 'SNPs'))
    all_or_some = st.radio("Select which subset of omics to use", ('All available omics', 'NDD-related omics', 'Custom'), help = 'NDD-related omics consist of omics related to brain areas,blood,and nerves.')
    remove_HLA = st.radio("Filter out all HLA genes?", ('Remove HLA genes (default)', 'Keep HLA genes'))

    if remove_HLA == 'Remove HLA genes (default)':
        HLA_status = True
    else:
        HLA_status = False
        remove


    submitted = st.form_submit_button("Next")
    if submitted:
        st.session_state['query_submit'] = 'run'

if st.session_state['query_submit'] == 'run':
    with st.form('parameter_selection'):
        # disease choice
        unique_dx = list(main_df['Disease'].unique())
        unique_dx.append('All')

        diseases = st.multiselect('Please select disease(s)', unique_dx)

        

        if "All" in diseases:
            diseases = list(main_df['Disease'].unique())
            st.session_state['dx_list'] =  diseases
            #st.write(f'**Selected Diseases**: {", ".join(st.session_state["dx_list"])}')
        else:
            st.session_state['dx_list'] =  diseases
            #st.write(f'**Selected Diseases**: {", ".join(st.session_state["dx_list"])}')

        if all_or_some == 'Custom':
        # custom omics_list
            unique_omic = list(main_df['Omic'].unique())

            omics = st.multiselect('Please select omic(s)', unique_omic)
            st.session_state['omic_list']  = omics
        

        # parameter choice
        pvalue = st.number_input('SMR P-value threshold', value = 0.05)
        st.session_state['pval'] = pvalue
        heidival = st.number_input('HEIDI p-value threshold', value = 0.01)
        st.session_state['heidi'] = heidival
        peqtlval = st.number_input('p-eQTL threshold', value = 0.05)
        st.session_state['peqtl'] = peqtlval

        submitted = st.form_submit_button("Submit parameters")
        if submitted:
            st.session_state['top_submit'] = 'run'

        # run fuctions once user hits run!
        if st.session_state['top_submit'] == 'run':

            # User choses to extract top genes across all omics
            if gene_or_not == 'Genes' and all_or_some == 'All available omics':
                st.write(f'Calculating results for Genes + All available omics')
                    
                top_genes_df, top_gene_chart = get_top_genes(main_df, st.session_state['dx_list'],
                omics_list = 'all', HLA = HLA_status, create_chart = True,heidi = st.session_state['heidi'],
                pval = st.session_state['pval'], peqtl = st.session_state['peqtl'])
                
                st.session_state['top_df'] = top_genes_df
                st.session_state['top_chart'] = top_gene_chart

            # User choses to extract top genes across NDD-related omics  
            elif gene_or_not == 'Genes' and all_or_some == 'NDD-related omics':
                st.write('Calculating results for Genes + NDD-related omics')
                    
                top_genes_df, top_gene_chart = get_top_genes(main_df, st.session_state['dx_list'],
                HLA = HLA_status, create_chart = True, heidi = st.session_state['heidi'],
                pval = st.session_state['pval'], peqtl = st.session_state['peqtl'])

                st.session_state['top_df'] = top_genes_df
                st.session_state['top_chart'] = top_gene_chart
            
            # User choses to extract top genes across custom omics list
            elif gene_or_not == 'Genes' and all_or_some == 'Custom':
                st.write(f'Calculating results for Genes + {", ".join(st.session_state["omic_list"])}')

                top_genes_df, top_gene_chart =  get_top_genes(main_df, st.session_state['dx_list'],
                omics_list = omics, HLA = HLA_status, create_chart = True, heidi = st.session_state['heidi'], 
                pval = st.session_state['pval'], peqtl = st.session_state['peqtl'])

                st.session_state['top_df'] = top_genes_df
                st.session_state['top_chart'] = top_gene_chart

            # User choses to extract top SNPs across all omics
            elif gene_or_not == 'SNPs' and all_or_some == 'All available omics':
                st.write('Calculating results for SNPs + All available omics')
                
                top_snps_df, top_snps_chart = get_top_snps(main_df, st.session_state['dx_list'],
                    omics_list = 'all', HLA = HLA_status, create_chart = True, heidi = st.session_state['heidi'],
                    pval = st.session_state['pval'], peqtl = st.session_state['peqtl'])

                
                st.session_state['top_df'] = top_snps_df
                st.session_state['top_chart'] = top_snps_chart
        
            # User choses to extract top SNPs across NDD-related omics    
            elif gene_or_not == 'SNPs' and all_or_some == 'NDD-related omics':
                st.write('Calculating results for SNPs + NDD-related omics')
                
                top_snps_df, top_snps_chart = get_top_snps(main_df, st.session_state['dx_list'], 
                HLA = HLA_status, create_chart = True, heidi = st.session_state['heidi'], 
                pval = st.session_state['pval'], peqtl = st.session_state['peqtl'])

                st.session_state['top_df'] = top_snps_df
                st.session_state['top_chart'] = top_snps_chart

            # User choses to extract top SNPs across acustom omics list   
            elif gene_or_not == 'SNPs' and all_or_some == 'Custom':
                st.write(f'Calculating results for SNPs + {", ".join(st.session_state["omic_list"])}')
                
                top_snps_df, top_snps_chart = get_top_snps(main_df, st.session_state['dx_list'], 
                omics_list = omics, HLA = HLA_status, create_chart = True, heidi = st.session_state['heidi'], 
                pval = st.session_state['pval'], peqtl = st.session_state['peqtl'])

                st.session_state['top_df'] = top_snps_df
                st.session_state['top_chart'] = top_snps_chart

            else:
                st.write('Error: Not a valid option')


with st.container():
    if st.session_state['top_submit'] == 'run':
        st.header('Result Dataframe')
        st.dataframe(st.session_state['top_df'])

        # download 
        output_name = st.text_input('Please provide an output file name if you would like to download your results dataframe', placeholder = 'example.csv')
        st.session_state['df_name'] = output_name

        if st.session_state['df_name']:
            st.download_button(label="Download data as CSV", data=convert_df(st.session_state['top_df']),file_name=st.session_state['df_name'], mime='text/csv')

        st.header('Result Chart')
        st.dataframe(st.session_state['top_chart'])

        chart_name = st.text_input('Please provide an output file name if you would like to download your results chart', placeholder = 'example.csv')
        st.session_state['chart_name'] = chart_name
            
        if st.session_state['chart_name']:
            st.download_button(label="Download data as CSV", data=convert_df(st.session_state['top_chart']),file_name=st.session_state['chart_name'], mime='text/csv')

with st.container():
    col1, col2= st.columns(2)

    with col1:
        with st.form('druggable'): # form to compare users df to druggable genes list
            st.subheader("Druggable gene identification")
            st.write('Identify genes in your top genes chart that are known drug targets')
            
            drug_run = st.form_submit_button("Run!")
            
            
            if drug_run: # if run selection is pressed, run adjustment on filtered df

                # calculate value counts
                value_df = gene_count(st.session_state['top_chart'])
                st.session_state['value_df'] = value_df
                st.session_state['value_status'] = 'run'
                drug_genes = drug_identify(value_df,st.session_state['drugdf'])
                st.session_state['gene_drugdf'] = drug_genes
                st.session_state['drug_status'] = 'run'

        if st.session_state['drug_status'] == 'run':
            st.dataframe(st.session_state['gene_drugdf'])

            # list of druggable genes
            gene_list = list(st.session_state['gene_drugdf']['gene'].unique())
            st.write(f'The druggable genes in your data are: {", ".join(gene_list)}')

            drug_filename = st.text_input('Please provide an output file name for the druggable gene results', placeholder = 'example.csv')
            st.session_state['drug_name'] = drug_filename
            
            if st.session_state['drug_name']:
                st.download_button(label="Download data as CSV", data=convert_df(st.session_state['gene_drugdf']),file_name=st.session_state['drug_name'], mime='text/csv')
    
    with col2:
        with st.form('value counts'):
            st.subheader("Gene/SNP counts")
            st.write('This tool will return the counts of how many times a gene or SNP shows up in your top results.')

            value_run = st.form_submit_button("Return value counts")

            if value_run:
                if st.session_state['value_status'] == 'run':
                    st.dataframe(st.session_state['value_df'])
                    top = st.session_state['value_df'].iloc[0,0]
                    st.write(f'The top results is {top}')
            
                else:
                    value_df = gene_count(st.session_state['top_chart'])
                    st.session_state['value_df'] = value_df
                    st.session_state['value_status'] = 'run'

                    st.dataframe(st.session_state['value_df'])
                    top = st.session_state['value_df'].iloc[0,0]
                    st.write(f'The top results is {top}')
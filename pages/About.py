import streamlit as st

st.title('OmicSynth')

st.write('OmicSynth aims to provide users with the ability to navigate Genome Wide Assocation Study (GWAS) summary statistics and SMR Analysis for multiple neurodegenerative diseases.')
st.write("All data used is publicly available")
st.write("GWAS summary statistics used include:")
st.write("""
1. AD - [Bellenguez 2022]('https://www.nature.com/articles/s41588-022-01024-z')
2. ALS - [Nicolas 2019]('https://pubmed.ncbi.nlm.nih.gov/29566793/')
3. FTD - [Ferrari 2014]('https://pubmed.ncbi.nlm.nih.gov/24943344/')
4. LBD - [Chia 2021]('https://www.nature.com/articles/s41588-021-00785-3')
5. PD - [Nalls 2019]('https://pubmed.ncbi.nlm.nih.gov/31701892/')
6. PSP - [Hoglinger 2011]('https://pubmed.ncbi.nlm.nih.gov/21685912/') """)
st.write("Reference Panel: 1000 Genomes")
st.write("Omic Data Sources: GTEx v8, metaBrain, ... available in SMR ready format [here](https://yanglab.westlake.edu.cn/software/smr/#DataResource)")
st.write("""1. GTEx v8 eQTL data 
- [Original Source]('https://www.gtexportal.org/home/datasets') 
- [SMR formatted and HG19/GR37 source]('https://yanglab.westlake.edu.cn/software/smr/#DataResource')""")
st.write("Druggable genome list source: [Finan et al., 2017]('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6321762/')")
st.write("Data curated by NIH CARD and Data Tecnica International.")
st.write("Github with code: ")
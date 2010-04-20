DEBUG= True

#DATABASE Settings
DB_SERVER = "gobi.mssm.edu"
DB_USER = "reader"
DB_PASSWORD = "sseltoor"
DB_DATABASE = "geneius"

#GENOME 2bit files
GENOME_PATH="/bio/ucsc/%/%.2bit"

#NCBI WSDL Files
#if relative path, make it in 
#reference to the geneius.py script,
#that is who will be importing it 

EUTILS_WSDL="http://www.ncbi.nlm.nih.gov/entrez/eutils/soap/v2.0/eutils.wsdl"
SEQUENCES_WSDL="http://www.ncbi.nlm.nih.gov/entrez/eutils/soap/v2.0/efetch_seq.wsdl"

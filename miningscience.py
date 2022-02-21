#Vamos a crear la funciones download_pubmed y mining_pubs para poder descargar un dataframe de Pubmed
from Bio import Entrez 
from Bio import SeqIO
import pandas as pd
import re,csv,itertools
import numpy as np

#Se define la función de download_pubmed para tener un dataframe con los ID de los artículos con una keyword específica
def download_pubmed(keyword):
    """"La funcion download_pubmed está encargada de recuperar la ID de los artículos científicos que esten relacionados con la palabra clave que se envió en la base de datos Pubmed """
    Entrez.email = 'leonardo.proano@est.ikiam.edu.ec'
    handle = Entrez.esearch(db='pubmed',retmax=10**5 ,retmode='xml',term=keyword)
    results = Entrez.read(handle)
    return results

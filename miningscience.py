#Vamos a crear la función download_pubmed para poder descargar un dataframe de Pubmed
from Bio import Entrez 
from Bio import SeqIO
import pandas as pd
import re,csv,itertools
import numpy as np

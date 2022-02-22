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

#Ahora, definimos la función mining_pubs    
def mining_pubs(tipo):
    """"La función de mining_pubs tiene la tarea de filtrar la data dependiendo de la variable (tipo) como DP AU o AD. Por otro lado, los datos  PMID para la extracción de archivos se obtienen gracias a la función download_pubmed, donde ecuentra los documentos exactos basados en su ID   """
    # Uso del metodo download_pubmed
    results = download_pubmed('Ecuador genomics') 
    
    
    id_list = results['IdList']                                                  #separamos IDs    
    ids = ','.join(id_list)    
    Entrez.email = 'leonardo.proano@est.ikiam.edu.ec'
    handle = Entrez.efetch(db='pubmed',rettype='medline',retmode='text',id=ids)  #rescatamos ducumentos por cada ID
    all_data = handle.read()                                                     #Revisamos la data
    
    if(tipo == "DP"):                                    #PMID y DP_year 
        zipcodes = re.findall(r'PMID-.(.+)', all_data)
        zipcodes1 = re.findall(r'DP  -.(.+)', all_data)
        all_ = list(zip(zipcodes,zipcodes1))
        nom_colum = ['PMID','DP_year']
    else:
        if(tipo == "AU"):                                #PMID y num_auth
            zipcodes = re.findall(r'PMID-.(.+)|(AU)  -|', all_data) 
            nom_colum = ['PMID','num_auth']
            
        elif(tipo == "AD"):                              #country y num_auth
            zipcodes = re.findall(r'PL  -.(.+)|(AU)  -|', all_data)
            nom_colum = ['country','num_auth']
        mira = list()
        for x in zipcodes:
            if(x[0]!=''):
                mira.append((x[0],''))
            elif(x[1]!=''):
                mira.append(('',x[1]))
        zipcodes= mira       
        lista_1 = list()
        lista_2 = list()
        va_c = 0
        for y in zipcodes:
            if(y[0] !=''):
                x_0 = y[0]
                lista_1.append(y[0])
                if(va_c != 0):
                    lista_2.append(va_c)
                    va_c = 0
            else:
                va_c = va_c+1            
        all_ = list(zip(lista_1,lista_2))
        
    results = pd.DataFrame(all_,columns = nom_colum)             
    return results

#Aquí se enccuentra la sección de ejecución del codigo y envio de variables
#id_list = results['IdList']          #separamos IDs 

if __name__ == '__main__':
 #--------------------------Se ingresa la variable tipo para iniciar---------------------------------------
    
    #resultado_final = mining_pubs("AU")  #Enviamos el tipo para el procesamiento ER (DP AU AD)  

#---------------------------------------------------------------------------------------------------------
  
    print("El nombre de la la funcion es:",download_pubmed.__name__ )
    print("Documentacion de la funcion :",download_pubmed.__doc__)
    print("____________________________________________________________________________________________________________________________")
    print("El nombre de la la funcion es:",mining_pubs.__name__ )
    print("Documentacion de la funcion :",mining_pubs.__doc__)


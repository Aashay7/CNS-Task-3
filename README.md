# CNS Task - 3 : Retrieve CellxGene Data

#### _by Aashay Gondalia_ (aagond@iu.edu)

In this workbook, the private download API of cellxgene website is used to retrieve a table for all "Normal" and "Homo Sapiens" data with following fields:
- Organ/Tissue Type
- Cell Type CL ID
- HGNC/ENSEMBL Gene IDs separated by ;
- No of. cells of this type
- Disease
- Assay
- Tissue
- Dataset Name

The table has been prepared in the pipe("|") seperated values format. 

Scanpy is used to read the h5 format data. With the help of this package, the data is accessed and used to populate the table seamlessly.

> NOTE: Due to RAM limitation, there were issues with accessing some of the larger datasets. 
Issue
> Scanpy was unable to read the dataset, which resulted in the vital information being inacessible for the table. 
> The work was executed on Google Colab and the RAM limitation of 12GB hindered the data reads on the larger datasets.
> NOTE: On the cellxgene website, there are some datasets with an empty cell_count field. 

> For results, I have attached the screenshots and images of the table(in the pipe-seperated value format as well as DataFrame format)

## Code

#### Importing necessary packages


  
```sh
import datetime
import requests
from requests.adapters import HTTPAdapter
import json
import numpy as np
import pandas as pd
import scanpy as sc
import os
```

#### Function to fetch the Collection Data from the https://cellxgene.cziscience.com/ website
```sh
def fetchCollectionData():
  ## HTTP Adapter Setup
  adapter = HTTPAdapter(max_retries=3)  #Hard-coded 3 Max Retries
  https = requests.Session()
  https.mount("https://", adapter)

  ## URL Elements
  CELLXGENE_PRODUCTION_ENDPOINT = 'https://api.cellxgene.cziscience.com'
  COLLECTIONS = CELLXGENE_PRODUCTION_ENDPOINT + "/dp/v1/collections/"
  DATASETS = CELLXGENE_PRODUCTION_ENDPOINT + "/dp/v1/datasets/"

  ## Fetch collection data
  r = https.get(COLLECTIONS)
  r.raise_for_status()

  collections = sorted(r.json()['collections'], key= lambda key :key['created_at'], reverse=True)
  print('Collection Fetch Complete.')
  return collections, https, CELLXGENE_PRODUCTION_ENDPOINT, COLLECTIONS, DATASETS
```


#### Function to apply filter while accessing the datasets - Applied Filters : {'Disease': 'Normal', 'Species': 'Homo Sapiens'}

```sh
def filter_Dataset_Homo_Sapien_Normal(all_collections):
  only_normal_homo_sapiens_ids = []
  for metadata in all_collections:
    collection_cell_counter = 0
    for dataset in metadata['datasets']:
      diseases = dataset['disease']
      id = dataset['id']
      for disease in diseases:
        if (str(disease['label']).lower() == 'normal' and str(dataset['organism']['label']).lower() == 'homo sapiens'):
          #Disease = disease['label']
          #Assay = dataset['assay']
          #Tissue = dataset['tissue']
          #Dataset_Name = dataset['name']
          try:
            collection_cell_counter += dataset['cell_count']
          except:
            pass
          only_normal_homo_sapiens_ids.append(id)
  print('Dataset Filters Applied.')
  return only_normal_homo_sapiens_ids
```

#### Initializing the output table. As mentioned in the google document, the dataframe column names are set accordingly. 

```
def initializeTable():
  table = pd.DataFrame({
    'Organ/Tissue Type' : [], 
    'Cell Type CL ID' : [],
    'HGNC/ENSEMBL Gene IDs' : [],
    'No. of Cells of this type' : [],
    'Disease' : [],
    'Assay' : [],
    'Tissue' : [],
    'Dataset Name' : [],
    })
  print('Table Initialization Complete.')
  return table
```

#### Gene data is fetched from the X(embeddings) matrix and all the non-zero values of GENEs for the corresponding cells are fetched.

```
def fetch_and_include_GENE_data(table, dataset, list_ex):
  print('Introducting GENE data into the table.')
  rows,cols = dataset.X.nonzero()
  #print('Rows - ', rows , len(rows))
  #print('Cols - ', cols, len(cols))

  diction = {}
  for i in list_ex:
    diction[i] = []

  #print(diction)

  for j in range(len(cols)):
    if (rows[j] in diction.keys()):
      diction[rows[j]].append(cols[j])

  gene_list = []
  #dataset_rows = len(diction.keys())
  for i in diction.keys():
    genes = ''
    for position in diction[i]:
      genes = genes + dataset.var_names[position] + ';'
    gene_list.append(genes)
  #print(len(gene_list))
  #print(gene_list)
  #for g in gene_list:
  #  print(g, '\n')
  for i in range(table.shape[0]):
    table['HGNC/ENSEMBL Gene IDs'][i] = gene_list[i]
  print('Gene Data Succesfully entered.')
```

#### This function is used to read the downloaded data from the cellxgene website, fetch the required information and enter it into the table. 

```
def enter_Details_into_Table(download_name, Disease, Assay, Tissue, Dataset_Name):
  '''try:
    table = pd.read_csv('dataTable.csv', sep='|')
    print('Table already exists -> Imported Data')
  except:'''
  table = initializeTable()

  print('Adding data to Table....')
  dataset = sc.read_h5ad(download_name)
  print('Dataset Imported in Scanpy Successfully.')
  os.remove(download_name)
  print('Removed dataset to aid program execution.')

  #print('Dataset Reading Complete')
  # 'Organ/Tissue Type', 'Cell Type CL ID', 'HGNC/ENSEMBL Gene IDs',
  # 'Cells of this type', 'Disease', 'Assay', 'Tissue', 'Dataset Name'
  
  # Gene IDs Aggregation into a single field.
  table_cell_ids = []
  table_row_ids = []
  for i in range(dataset.shape[0]):

    if (dataset.obs['cell_type_ontology_term_id'][i] not in table_cell_ids):
      table_row_ids.append(i)
      
      num_cells_float = dataset.obs.cell_type_ontology_term_id.value_counts()[dataset.obs['cell_type_ontology_term_id'][i]] 
      no_of_cells_of_same_type = int (num_cells_float)
      
      table.loc[len(table.index)] = [
                              dataset.obs['tissue'][i],
                              dataset.obs['cell_type_ontology_term_id'][i],
                              '_to_be_filled_',
                              no_of_cells_of_same_type,
                              Disease, 
                              dataset.obs['assay'][i], 
                              Tissue, 
                              Dataset_Name
                              ]
      
      table_cell_ids.append(dataset.obs['cell_type_ontology_term_id'][i])

  fetch_and_include_GENE_data(table, dataset, table_row_ids)
  print('Data successfully added to the Table.')
  print(table)

  return table
```

#### Total checker function to match the total cell_count in the dataset and the cell_count mentioned on the website.

```
def check_total_cell_count(cell_count_total_website, table):
  print('Website Cell Count : ', cell_count_total_website)
  table_total = int(sum(table['No. of Cells of this type']))
  print('Total Cell Counts in the Table : ', table_total)
  try:
    if (table_total == cell_count_total_website):
      print('Cell count in the prepared table matches the cell count data on the website!!')
  except:
    pass

```

#### Master Function is the main executable function. It calls all the above mentioned functions and saves the table in the required 'pipe-seperated' values format.

```
table_dataholder = None
def masterFunction():
  
  collections, https, CELLXGENE_PRODUCTION_ENDPOINT, COLLECTIONS, DATASETS = fetchCollectionData()
  all_collections = []

  ## INITIAL METADATA FETCH
  for collection in collections:
    r1 = https.get(COLLECTIONS + collection['id'], timeout=5)
    collection_metadata = r1.json()
    all_collections.append(collection_metadata)
  
  ## Populating only_normal_homo_sapiens_ids list with all the filtered dataset ids.
  only_normal_homo_sapiens_ids = filter_Dataset_Homo_Sapien_Normal(all_collections)
          
  
  for collection in all_collections:
    for dataset in collection['datasets']:
      try:
        cell_count_total_website = dataset['cell_count']
      except:
        cell_count_total_website = None

      for asset in dataset['dataset_assets']:

        # Using the H5 format for less overload and compatibility
        # Faced some issues with the RDS formatting. 
        #High overload on the python wrapper to read RDS files.
        if ((asset['filetype'] == 'H5AD') and (asset['dataset_id'] in only_normal_homo_sapiens_ids)):
          DATASET_REQUEST = DATASETS + asset['dataset_id']  +"/asset/"+  asset['id']
          
          r2 = requests.post(DATASET_REQUEST)
          r2.raise_for_status()
          presigned_url = r2.json()['presigned_url']
          
          headers = {'range': 'bytes=0-0'}
          r3 = https.get(presigned_url, headers=headers)
          print('\nDataset -> ', dataset['name'], '\nURL -> ', presigned_url)
          
          if (r3.status_code == requests.codes.partial):
            download_name = dataset['name'] + '.h5ad'
            print('Dataset Download Started.')
            r3 = https.get(presigned_url, timeout=10)
            r3.raise_for_status()
            open(download_name, 'wb').write(r3.content)
            print('Dataset Download Complete.')
            table = enter_Details_into_Table(download_name, dataset['disease'], dataset['assay'], dataset['tissue'], dataset['name'])

            ## SAVE TABLE 
            print('Saving Table....')
            os.mkdir(collection['name'])
            filepath = collection['name'] + '/' + dataset['name']
            table.to_csv(filepath, sep='|')
            print('Table Saved Successfully.')
            table_dataholder = table
            check_total_cell_count(cell_count_total_website, table)

            ## To effectively use the Google COLAB RAM. 
            table = None
            del table
            print('Local Copy of table removed from RAM')


```

#### Master Function currently fetches collection data and downloads the dataset in a serial manner, which can be parallelized.

```
masterFunction()
```

> For complete cell-by-cell execution results, I have also attached the ipynb file that was used for execution.
> The screenshots have been attached in this folder as well.

> The table size for a single dataset in pipe-seperated format is in GBs. The GeneID field is memory heavy field and it increases the overall size of the output table dataset.

### Screenshots 

![alt text](https://github.com/Aashay7/CNS-Task-3/blob/main/Screenshots/Capture.PNG)

#

![alt text](https://github.com/Aashay7/CNS-Task-3/blob/main/Screenshots/Screenshot%202021-08-22%20at%201.51.29%20PM.png)

# 

![alt text](https://github.com/Aashay7/CNS-Task-3/blob/main/Screenshots/Screenshot%202021-08-22%20at%203.26.52%20PM.png)


# 

> File Size Issue - Table size too big because of the Gene ID data field. Unable to attach the whole file in the repository.

![alt text](https://github.com/Aashay7/CNS-Task-3/blob/main/Screenshots/filesizeissue.PNG)



# Odor Meso

+ DataJoint pipeline to capture the data for the mice two-photon imaging experiments in the presence of olfactory stimulation.

+ This is a cloned of the [cajal/pipeline](https://github.com/cajal/pipeline) with modifications. 

---
## Install pipeline

### Step 1 - Clone the repositories

+ Launch a new terminal and run the following
    ```
    cd ~/<Enter custom path>
    git clone https://github.com/vathes/odor_meso
    git clone https://github.com/vathes/DataPlot
    cd odor_meso
    ```

### Step 2 - Create `.env`

+ At the root of the `odor_meso` repository folder, create a new file `.env` from the `env_template`

+ Specify the database's `DJ_HOST`, `DJ_USER`, and `DJ_PASS`

### Step 3 - Build the Docker image

+ Install [Docker Desktop](https://www.docker.com/products/docker-desktop)

+ From the root of the cloned `odor_meso` repository directory
     ```
     docker-compose up -d --build notebook
     ```

### Installation complete

+ At this point the basic setup of this pipeline is complete.
---
## Docker commands
+ The following is a resource for docker commands you may need to utilize.

+ Build the Docker image
     ```
     docker-compose up -d --build notebook
     ```

+ If the image is already built, start a container
     ```
     docker-compose up -d notebook
     ```

+ Access an interactive bash shell of the container
     ```
     docker exec -it odor_meso_notebook_1 /bin/bash
     ```

+  Stop and remove containers and networks
     ```
     docker-compose down
     ```
---
## Ingest data

+ Sample data must be ingested if you are using a local database

+ The database prefix is `pipeline`.  If you are running locally, ensure that you are not already using this database prefix.

+ For the DataMan function, place a `paths.init` file  in the following directory
     ```
     ~/odor_meso/external/dataplot_storage/
     ```

+ Place the fetched files in the following directory
     ```
     ~/odor_meso/external/ingestion_storage/animal_id_<animal_id>
     ```

+ The ingestion routine can be found within the following jupyter notebook.  To run the notebook, navigate to `localhost:8888` in your web broswer.
     ```
     ~/odor_meso/python/scripts/data_ingest.ipynb
     ```

+ Rerun the jupyter notebook when ingesting new scans and/or subjects

---
## Explore data

+ Detailed instructions can be found within the following jupyter notebook.  To run the notebook, navigate to `localhost:8888` in your web broswer.
     ```
     ~/odor_meso/python/scripts/data_explore.ipynb
     ```

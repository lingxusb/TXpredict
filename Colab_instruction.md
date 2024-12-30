## Instructions for using the Colab notebook to predcit transcriptome

### 1. connect to a GPU instance
click the triangle on the top right corner and select **Change runtime type**

<img width="290" alt="image" src="https://github.com/user-attachments/assets/66a95466-0667-463c-813a-d288b3762b68" />


Then select T4 instance and click **save**

<img width="454" alt="image" src="https://github.com/user-attachments/assets/6bcac394-5713-4c7e-bc95-2fd0a3e3b94a" />


After connecting to a T4 instance, you will see this on the top right corner:

![image](https://github.com/user-attachments/assets/3f7d8a53-e95c-4cf4-ab56-93b64a5149fd)

### 2. upload genome sequence and annotation files
click the **files** button on the left:

<img width="94" alt="image" src="https://github.com/user-attachments/assets/9c418843-5b75-4ef9-a274-3a0c5e858b62" />

Then click the **upload to session storage** button:

<img width="283" alt="image" src="https://github.com/user-attachments/assets/3d55c5a0-4a45-43ad-ae26-2c154971e2e7" />

Select genome sequence and annotation files from the local folders, please note that two files should have the same name: e.g., GCF_000011065.1.fna and GCF_000011065.1.gff

### 3. running the codes
Please click the **run** button for three parts of the codes, including **Install**, **Generate embeddings** and **Predict gene expression**:

<img width="225" alt="image" src="https://github.com/user-attachments/assets/0c395a7d-3359-4590-bd45-9797851fa26c" />

The predicted gene expression will be automatically downloaded as a csv file.


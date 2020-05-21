import os
import subprocess
import sys
import pandas as pd
import logging
import csv
import shutil
import smtplib
import mimetypes
from email.message import EmailMessage
import warnings


script=sys.argv[0]
base_dir=sys.argv[1]
database_dir=sys.argv[2]

def RunBlast(): 
     warnings.filterwarnings("ignore")
     os.chdir(base_dir)
     child_processes=[]
     os.system("mkdir blastoutput")
     for query in os.listdir('.'):
         if (query.endswith(".fasta") or (query.endswith(".fna"))):
             query=os.path.join(base_dir,query)
             queryfile=open(query)
             baseq=os.path.basename(query)
             filename =os.path.splitext(baseq)[0] 
             for database in os.listdir(database_dir):
                 database=os.path.join(database_dir,database)
                 basedb=os.path.basename(database)
                 dbname=basedb.split(".")[0]
                 databasename =os.path.join(database_dir,basedb.split(".")[0])
                 p=subprocess.Popen(["blastn","-query",query,"-db",databasename,"-evalue","1e-6","-outfmt","6 qseqid sseqid pident qlen slen qstart qend sstart send","-max_target_seqs","3","-out","./blastoutput"+"/"+filename+"_"+dbname+".blast"])
                 child_processes.append(p)
                 for cp in child_processes:
                     cp.wait()
#RunBlast()
#print("blast is done")


def filter():
     os.chdir(base_dir+"/blastoutput")
     os.system("mkdir sorted_blast_pair") 
     for blastresult in os.listdir('.'):
         if blastresult.endswith(".blast"):
             genomename=os.path.basename(blastresult)
             genomename=genomename.split(".")[0]
             blastresult=open(blastresult)
             for line in blastresult:
                 try:
                     gene={}
                     line = line.split( )
                     qseqid=line[0]
                     sseqid=line[1]
                     pident=float(line[2])
                     qlength=float(line[3])
                     slength=float(line[4])
                     qstart=float(line[5])
                     qend=float(line[6])
                     sstart=float(line[7])
                     send=float(line[8])
                    
                     if (pident> 97) & (abs(qend-qstart)/slength > 0.75)  :
                         gene[qseqid]=sseqid
                         for key in gene:
                             with open("./sorted_blast_pair" +"/"+genomename+".blast","a") as ofile:
                                  ofile.write(genomename+"\t"+key+"\t"+gene.get(key)+"\t"\
                                      +str(pident)+"\t"+str(slength)+"\t"+str(abs(qend-qstart)/slength)+"\t"+str(qstart)+"\t"+str(qend)+"\n")
                                  ofile.close 
                 except IOError:
                     print("no input")
             blastresult.close() 
#filter()
#print("Filtering blast result is done")


#### Generate tables for each genome"
def generate_table():
    os.chdir(base_dir+"/blastoutput/sorted_blast_pair")
    if not "result_tables" in os.listdir("."): 
            os.system("mkdir result_tables")
    for filename in os.listdir("."):
         if filename.endswith(".blast"):
             base=filename.split(".")
             name=base[0]
             df=pd.read_csv(filename,sep="\t",header=None)
             df.columns=['Genome','testgenome_gene_name','db_gene_name','pident','slength','coverage','querystart','queryend']
             df['db']=df.db_gene_name.apply(lambda x:x.split("_")[-1])
             df['score']=df['pident']*df['coverage']
             #df['result']=df.sort_values(by="score",ascending=False).head(1).db_gene_name.values[0]+"_true"
             #df=df.sort_values(by="score",ascending=False).head(6)
             f_name=base_dir+"/blastoutput/sorted_blast_pair/result_tables/"+str(filename)+"_table.csv"
             df.to_csv(f_name)

def result_table1():
     os.chdir(base_dir+"/blastoutput/sorted_blast_pair/result_tables")
     filelist=glob.glob("./*_table.csv")
     df_list=[pd.read_csv(file) for file in filelist]
     bigdf1=pd.concat(df_list,axis=0)
     bigdf1=bigdf1.drop("Unnamed: 0",axis=1)
    #bigdf=bigdf.drop("Unnamed: 0",axis=1)
    #  bigdf.index=['18s','actin','hsp70']
    #  bigdf.loc['18s']=bigdf.loc['18s'].str.replace("C.","C_")
    #  bigdf.loc['result']=''
     bigdf1.to_csv('result_table1.csv')

##filter best hit and insert newtype for not in database
def filter2():
    os.chdir(base_dir+"/blastoutput/sorted_blast_pair/result_tables")
    for file in os.listdir("."):
        if file.endswith("_table.csv"):
            dic1={}
            #dic2={}
            dic={}
            #dic_sub={}
            #result2=pd.DataFrame()
            dic1[file]=dic
            #dic2[file]=dic_sub
            sample=pd.read_csv(file)
            sample=sample.drop("Unnamed: 0",axis=1)
            sample['score']=sample['pident']*sample['coverage']
            table=pd.DataFrame(columns=sample.columns)
    #a.setdefault(sample.Genome, [])
            for i in ['18s','actin','hsp70']:
                    if i in str(sample['db']):
                        if (sample['db'].str.count(i).sum()) == 1:
                            dic[i]=sample[sample.db==i].db_gene_name.values[0]
                            table=table.append(pd.DataFrame(sample[sample.db_gene_name==str([dic[i]][0])]))
                        elif (sample['db'].str.count(i).sum()) > 1:
                            dic[i]=sample[sample.db==i].sort_values(by="score",ascending=False).head(1).db_gene_name.values[0]
                            table=table.append(pd.DataFrame(sample[sample.db_gene_name==str([dic[i]][0])].iloc[0]).T)
                            #print(dic_sub[i])
                    elif i not in str(sample['db']):
                        dic[i]="new type"
                        table=table.append(pd.DataFrame(sample[sample.db_gene_name==str([dic[i]][0])])).fillna("N/A")
            f_name="./"+str(file)+"_newtable.csv"
            f_name2="./"+str(file)+"_newtable2.csv"
            pd.DataFrame.from_dict(dic1).to_csv(f_name, quoting=csv.QUOTE_NONE,quotechar='',  escapechar=",")
            table.to_csv(f_name2)


def combine_table():
     os.chdir(base_dir+"/blastoutput/sorted_blast_pair/result_tables")
     import glob
     filelist=glob.glob("./*newtable.csv")
     df_list=[pd.read_csv(file) for file in filelist]
     bigdf=pd.concat(df_list,axis=1)
     bigdf=bigdf.drop("Unnamed: 0",axis=1)
     bigdf.index=['18s','actin','hsp70']
     bigdf.loc['result']=''
     for i in range(0,len(bigdf.columns)):
    #bigdf.loc['result']=bigdf[bigdf.columns[i]].apply(lambda x:"_".join(x.split("_")[0:2]).lower()).value_counts().index[0]
         bigdf.loc['result'][i]=bigdf[bigdf.columns[i]].apply(lambda x:"_".join(x.split("_")[0:2]).lower()).value_counts().index[0]
     bigdf.to_csv("species_call.csv")

     bigdf_T=bigdf.T.reset_index()
     bigdf_T.columns=['Genome','18s','actin','hsp70','result']
     bigdf_T.Genome=bigdf_T.Genome.apply(lambda x:x.split(".")[0])
     
     filelist1=glob.glob("./*_table.csv")
     filelist2=glob.glob("./*_newtable2.csv")

     df_list=[pd.read_csv(file) for file in filelist1]
     df_list2=[pd.read_csv(file) for file in filelist2]

     bigdf1=pd.concat(df_list,axis=0)
     bigdf2=pd.concat(df_list2,axis=0)
    #bigdf1=bigdf1.drop("Unnamed: 0",axis=0)
     bigdf2=pd.DataFrame(bigdf2.drop("Unnamed: 0",axis=1))
     bigdf1.to_csv("blast_merge.csv")
     bigdf2.to_csv("blast_merge2.csv")
     final=pd.merge(bigdf2,bigdf_T,on="Genome")
     final.to_csv("final_results.csv")


def outputdirectory():
     base_dir=sys.argv[1]
     os.chdir(base_dir)
     os.system("mkdir outputfiles")
     source_dir1=base_dir+"/blastoutput/sorted_blast_pair/result_tables"
     destination=base_dir+"/outputfiles"

     for file in os.listdir(source_dir1):
         if file.startswith("species_call"):
             target=os.path.join(destination,file)
             file=os.path.join(source_dir1,file)
             shutil.copy(file,target)
     for file in os.listdir(source_dir1):
         if file.startswith("final_results"):
             target=os.path.join(destination,file)
             file=os.path.join(source_dir1,file)
             shutil.copy(file,target)
     for file in os.listdir(source_dir1):
         if file.startswith("blast_merge"):
             target=os.path.join(destination,file)
             file=os.path.join(source_dir1,file)
             shutil.copy(file,target)

def sendemail():
     base_dir=sys.argv[1]
     os.chdir(base_dir+"/outputfiles")

     msg=EmailMessage()
     msg['Subject']="LpsubP result files"
     msg['From'] = "nej1@cdc.gov"
     msg['To'] = "nej1@cdc.gov",sys.argv[3]
     msg.preamble="Genotyper result files for"+str(sys.argv[1].split("/")[::-1][0])

     message = """
     Attached are the outputfiles from Genoyper.
     The general explanation for each file is also attached("output_explanation.txt").
     Let me know if you have any questions.
     Thanks,
     Subin
     """

     msg.set_content(message)

     for file in os.listdir("."):
         path=os.path.join(".",file)
         if not os.path.isfile(path):
             continue
         ctype,encoding = mimetypes.guess_type(path)
         if ctype is None or encoding is not None:
             ctype="application/octet-stream"
         maintype,subtype=ctype.split('/',1)
         with open(path,'rb') as fp:
             msg.add_attachment(fp.read(),maintype=maintype,subtype=subtype,filename=file)
            #msg.attach(fp.read())
        
         with smtplib.SMTP('localhost')as s:
              s.send_message(msg)


def main():
     logging.basicConfig(filename='myscripts.log',format='%(asctime)s %(message)s', level=logging.INFO)
     logging.info("Started")
     warnings.filterwarnings("ignore")
     RunBlast()
     warnings.filterwarnings("ignore")
     logging.basicConfig(format='%(asctime)s %(message)s')
     logging.info("Blast is done")
     filter()
     logging.basicConfig(format='%(asctime)s %(message)s')
     logging.info("First filtering is done")
     generate_table()
     filter2()
     logging.info("Second filtering is done")
     combine_table()
     logging.info("Result tables and final table were generated")
     outputdirectory()
     logging.info("outputdirectory is generated")
     sendemail()
     logging.info("results are sent by emails")
     logging.info('Finished')
#      code13()
#      logging.info('Outputfiles are sent to your email account')
if __name__ == '__main__':
     main()

# ####GetSequence#####

# os.chdir(base_dir)
# os.system("mkdir seqrecords")
# def Parse(filename,seqs):
#     file = open(filename) 
#     seqs={}
#     name = ''
#     for line in file:
#         line = line.rstrip()  
#         if line.startswith('>'):
#             name=line.replace('>',"")
#             seqs[name] = ''
#         else:
#             seqs[name] = seqs[name] + line
            
#     file.close
#     return seqs

# seqs={}    
# for genome in os.listdir('.'):

#      if genome.endswith(".fasta"):
#      #if genome.endswith(".fasta"):
#          seqs=dict(seqs,**Parse(genome,seqs))
            
# for file in os.listdir(base_dir+'/blastoutput/sorted_blast_pair'):
#      #genomename=file.split("_")[2]+file.split("_")[3]
#      genomename=file.split("_")[2]
#      #genomename=file.split(".")[0]
#      #genomename=genomename.split("_")[2]
#      file=open(os.path.join(base_dir+'/blastoutput/sorted_blast_pair',file))
#      for line in file:
#          genename=line.split("\t")[1]+" "
#          coregenename=line.split("\t")[0]
#          for key in seqs:   
#              if key.find(str(genename))!= -1: 
#              #if key == str(genename):
#                  with open("./seqrecords"+"/"+coregenename+"_"+genename+"_"+genomename+".fasta","w") as ofile:
#                       ofile.write(">"+coregenename+"_"+genename+"_"+genomename+"\n"+seqs.get(key))
#                       ofile.close()
#      file.close()
    
# print("Getting sequences are done")     

# os.chdir(base_dir+'/seqrecords')
# #os.chdir(base_dir)
# os.system('mkdir pergene_seqrecords')
# genelist=open(os.path.join(sys.argv[1],'new_49gene.list'))

# for gene in genelist: 
#      gene=gene.rstrip()
#      for seqrecord in os.listdir("."):
#          if seqrecord.startswith(gene):
#              seq=open(os.path.join(base_dir+'/seqrecords',seqrecord))
#              #seq=open(os.path.join(base_dir,seqrecord))
#              for seqline in seq:
#                  seqline=seqline.rstrip()
#                  #with open("./pergene_seqrecords/"+"/"+gene+"blast_list"+".txt","a") as pfile:
#                  #    pfile.write(seqline+"\n")
#                  #    pfile.close()
#                  with open("./pergene_seqrecords"+"/"+gene+"_"+"unaligned"+".fasta","a") as pfile:
#                       pfile.write(seqline+"\n")
#                       pfile.close
#              seq.close()
# genelist.close()
# print("Sequences are sorted by each locus")

# #####PRESENCE/ABSENCE#####
# os.chdir(base_dir)
# filelist1=os.listdir(base_dir+"/blastoutput")
# filelist2=os.listdir(base_dir+"/blastoutput/sorted_blast_pair")
# sys.stdout=open('test','a')
# for beforefile in filelist1:
#      if beforefile.endswith(".blast"):
#          base=beforefile.split(".")[0]
#          coregenename=base.split("_")[0]
#          #genomename=base.split("_")[1]+"_"+base.split("_")[2]
#          genomename=base.split("_")[1]
     
     
#          if str(filelist2).find(str(beforefile))!= -1:
#              sys.stdout.write(coregenename+"\t"+genomename+"\t"+"yes"+"\n")      
#          else:
#              sys.stdout.write(coregenename+"\t"+genomename+"\t"+"no"+"\n")
         
# sys.stdout.close()


# test=open("test")
# no=open("notest",'a')
# for line in test: 
#      line=line.rstrip()
#      core=line.split()[0]
#      subject=line.split()[1]
#      if (line.startswith("lpg")) & (line.find("no")!= -1):
#          for blastresult in filelist1:
#              if blastresult.startswith(core+"_"+subject):
#                  f=open(os.path.join(base_dir+"/blastoutput",blastresult))
#                  no.write(line+"\t"+str(f.readlines()).replace("t","").replace("n","")+"\n")

# no.close()

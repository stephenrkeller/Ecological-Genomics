## Author: Stephen R. Keller  
## Ecological Genomics:   

### Overall Description of notebook      

My notes on coding and data analysis for the Spring 2020 Ecological Genomics course @ UVM


### Date started: 2020-Jan-22
### Date end:   (year-month-day)    

### Philosophy   
Science should be reproducible and one of the best ways to achieve this is by logging research activities in a notebook. Because science/biology has increasingly become computational, it is easier to document computational projects in an electronic form, which can be shared online through Github.    

### Helpful features of the notebook     

**It is absolutely critical for your future self and others to follow your work.**     

* The notebook is set up with a series of internal links from the table of contents.    
* All notebooks should have a table of contents which has the "Page", date, and title (information that allows the reader to understand your work).     
* Also, one of the perks of keeping all activities in a single document is that you can **search and find elements quickly**.     
* Lastly, you can share specific entries because of the three "#" automatically creates a link when the notebook renders on github.      


<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.  


### Table of contents for 60 entries (Format is *Page: Date(with year-month-day). Title*)        
* [Page 1: 2020-01-22](#id-section1). Intro to Github, Markdown, and UNIX command-line
* [Page 2:](#id-section2).
* [Page 3:](#id-section3).
* [Page 4:](#id-section4).
* [Page 5:](#id-section5).
* [Page 6:](#id-section6).
* [Page 7:](#id-section7).
* [Page 8:](#id-section8).
* [Page 9:](#id-section9).
* [Page 10:](#id-section10).
* [Page 11:](#id-section11).
* [Page 12:](#id-section12).
* [Page 13:](#id-section13).
* [Page 14:](#id-section14).
* [Page 15:](#id-section15).
* [Page 16:](#id-section16).
* [Page 17:](#id-section17).
* [Page 18:](#id-section18).
* [Page 19:](#id-section19).
* [Page 20:](#id-section20).
* [Page 21:](#id-section21).
* [Page 22:](#id-section22).
* [Page 23:](#id-section23).
* [Page 24:](#id-section24).
* [Page 25:](#id-section25).
* [Page 26:](#id-section26).
* [Page 27:](#id-section27).
* [Page 28:](#id-section28).
* [Page 29:](#id-section29).
* [Page 30:](#id-section30).
* [Page 31:](#id-section31).
* [Page 32:](#id-section32).
* [Page 33:](#id-section33).
* [Page 34:](#id-section34).
* [Page 35:](#id-section35).
* [Page 36:](#id-section36).
* [Page 37:](#id-section37).
* [Page 38:](#id-section38).
* [Page 39:](#id-section39).
* [Page 40:](#id-section40).
* [Page 41:](#id-section41).
* [Page 42:](#id-section42).
* [Page 43:](#id-section43).
* [Page 44:](#id-section44).
* [Page 45:](#id-section45).
* [Page 46:](#id-section46).
* [Page 47:](#id-section47).
* [Page 48:](#id-section48).
* [Page 49:](#id-section49).
* [Page 50:](#id-section50).
* [Page 51:](#id-section51).
* [Page 52:](#id-section52).
* [Page 53:](#id-section53).
* [Page 54:](#id-section54).
* [Page 55:](#id-section55).
* [Page 56:](#id-section56).
* [Page 57:](#id-section57).
* [Page 58:](#id-section58).
* [Page 59:](#id-section59).
* [Page 60:](#id-section60).

------
<div id='id-section1'/>
### Page 1: 2018-01-24. Notes on using Github, markdown, and the UNIX command-line

* Ask Lauren to give a brief intro on using git to create a repo and document your work in a lab notebook; push to origin on github

* I'll lead a tutorial on logging into the class unix server, doing some basic unix navigation and file manipulation, and writing simple for loops to run Fastqc.  Notes for all the 2020 tutorials are posted here: [2020 tutorial page](https://pespenilab.github.io/Ecological-Genomics/Tutorials.html)

#### Let's get started!

*  Embed a bit of unix code showing how to log-in to the server:

```
ip0af52994:~ stephenkeller$ ssh srkeller@pbio381.uvm.edu
srkeller@pbio381.uvm.edu's password: 
Last login: Thu Jan 16 10:58:54 2020 from ip0af52994.int.uvm.edu
[srkeller@pbio381 ~]$ pwd
/users/s/r/srkeller
[srkeller@pbio381 ~]$ hostname
pbio381.uvm.edu
[srkeller@pbio381 ~]$ 
```

* Show where the project data  fastq files live:

```
[srkeller@pbio381 myresults]$ cd /data/project_data/RS_ExomeSeq/fastq/edge_fastq/
[srkeller@pbio381 edge_fastq]$ ls
AB_05_R1_fastq.gz   CRA_02_R2_fastq.gz  GPA_03_R1_fastq.gz  MRC_16_R2_fastq.gz  RP_11_R1_fastq.gz   XDS_06_R2_fastq.gz

[...]
```

How to list out just some of the pops (say, just the AB's), or create a list of the number of unique pops?

```
[srkeller@pbio381 edge_fastq]$ ls AB*
AB_05_R1_fastq.gz  AB_08_R1_fastq.gz  AB_12_R1_fastq.gz  AB_16_R1_fastq.gz  AB_18_R1_fastq.gz
AB_05_R2_fastq.gz  AB_08_R2_fastq.gz  AB_12_R2_fastq.gz  AB_16_R2_fastq.gz  AB_18_R2_fastq.gz

[srkeller@pbio381 edge_fastq]$ ls | sed 's/\_/\t/g' | cut -f1 | uniq >~/myresults/metadata/Edge.pops 
[srkeller@pbio381 edge_fastq]$ cat ~/myresults/metadata/Edge.pops 
AB
BFA
BRB
CR
CRA
CRR
DG
GFM
GPA
HR
KOS
MRC
MT
PRK
RP
WA
XCS
XCV
XDS
XFS
XGL
XPK
XSK
XWS
```

#### We can use `scp` or a GUI sftp app to transfer the files to our local drive and then upload to github. 

```
```

* Here's a picture from one of the FASTQC run to demonstrate how to embed links to pictures from your online lab notebook. ![](https://github.com/stephenrkeller/Ecological_Genomics/blob/master/Spring_2020/myresults/fastqc/AB_05_R1_fastq_fastqc/Images/per_base_quality.png) 



------
<div id='id-section2'/>
### Page 2: 2018-01-29 
------
<div id='id-section3'/>
### Page 3:
------
<div id='id-section4'/>
### Page 4:
------
<div id='id-section5'/>
### Page 5:
------
<div id='id-section6'/>
### Page 6:
------
<div id='id-section7'/>
### Page 7:
------
<div id='id-section8'/>
### Page 8:
------
<div id='id-section9'/>
### Page 9:
------
<div id='id-section10'/>
### Page 10:
------
<div id='id-section11'/>
### Page 11:
------
<div id='id-section12'/>
### Page 12:
------
<div id='id-section13'/>
### Page 13:
------
<div id='id-section14'/>
### Page 14:
------
<div id='id-section15'/>
### Page 15:
------
<div id='id-section16'/>
### Page 16:
------
<div id='id-section17'/>
### Page 17:
------
<div id='id-section18'/>
### Page 18:
------
<div id='id-section19'/>
### Page 19:
------
<div id='id-section20'/>
### Page 20:
------
<div id='id-section21'/>
### Page 21:
------
<div id='id-section22'/>
### Page 22:
------
<div id='id-section23'/>
### Page 23:
------
<div id='id-section24'/>
### Page 24:
------
<div id='id-section25'/>
### Page 25:
------
<div id='id-section26'/>
### Page 26:
------
<div id='id-section27'/>
### Page 27:
------
<div id='id-section28'/>
### Page 28:
------
<div id='id-section29'/>
### Page 29:
------
<div id='id-section30'/>
### Page 30:
------
<div id='id-section31'/>
### Page 31:
------
<div id='id-section32'/>
### Page 32:
------
<div id='id-section33'/>
### Page 33:
------
<div id='id-section34'/>
### Page 34:
------
<div id='id-section35'/>
### Page 35:
------
<div id='id-section36'/>
### Page 36:
------
<div id='id-section37'/>
### Page 37:
------
<div id='id-section38'/>
### Page 38:
------
<div id='id-section39'/>
### Page 39:
------
<div id='id-section40'/>
### Page 40:
------
<div id='id-section41'/>
### Page 41:
------
<div id='id-section42'/>
### Page 42:
------
<div id='id-section43'/>
### Page 43:
------
<div id='id-section44'/>
### Page 44:
------
<div id='id-section45'/>
### Page 45:
------
<div id='id-section46'/>
### Page 46:
------
<div id='id-section47'/>
### Page 47:
------
<div id='id-section48'/>
### Page 48:
------
<div id='id-section49'/>
### Page 49:
------
<div id='id-section50'/>
### Page 50:
------
<div id='id-section51'/>
### Page 51:
------
<div id='id-section52'/>
### Page 52:
------
<div id='id-section53'/>
### Page 53:
------
<div id='id-section54'/>
### Page 54:
------
<div id='id-section55'/>
### Page 55:
------
<div id='id-section56'/>
### Page 56:
------
<div id='id-section57'/>
### Page 57:
------
<div id='id-section58'/>
### Page 58:
------
<div id='id-section59'/>
### Page 59:
------
<div id='id-section60'/>
### Page 60:

------

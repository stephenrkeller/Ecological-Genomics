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
* [Page 2:](#id-section2). FastQC, read trimming, and mapping to reference genome
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



* Went through  a bit of unix code showing how to log-in to the server

* Showed where the shared project data files live

* Went through unix commands:  `pwd, cd, ls, grep, cp, wc, pipe (|), redirect to file (>)`

* Showed how to use git to add, commit, and push to server

* Posted to Slack afterwards about using `git pull` prior to merge any changes from github first prior to pushing commits to the cloud

* On next Wednesday, pick up with using vim to edit files (using .bashrc as example) and then working with fastq files

* Plan should be to walk-through fastqs, assign population samples to students, have them run fastqc, have them design a trimomatic script to clean the reads up, then set up a bash script in `screen` to map reads with bwa and calculate genotype likelihoods 


------
<div id='id-section2'/>
### Page 2: 2018-01-29 

Example output from FastQC (copied to my 'docs' folder on GitHub for webpage display)

![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAyAAAAJYCAIAAAAVFBUnAAAjV0lEQVR42u3dS3LbSBZAUex/M1qEJ95BT3pc49oBW2GF2Qh8Ml9+kACI80JR4ZJlHYoJgVcACU0vY4wxxhjTdSZ3gTHGGGOMwDLGGGOMEVjGGGOMMQLLGGOMMcYILGOMMcYYgWWMMcYYI7CMMcYYY4zAMsYYY4wRWMYYY4wxAssYY4wxxggsY4wxxhiBZYwxxhgjsIwxxhhjBJYxxhhjjBFYxhhjjDECyxhjjDFGYBljjDHGGIFljDHGGCOwjDHGGGMEljHX2HAnm+7ou3rzPj99IaY/86l3rxmzN7j7VmQEljH/35c17s7WD0v2j48NrM2t61Lbw/rGHLEBd/mSg/fepe7nLvehHYgRWOZzftbsuze0fxRY69h6VGAdeoPXdXWdm92rLH2zG4FlPqGu0o8ui0fN9c/Ki3+4+JjIA+3mz9/pd65vwOKRJv45IwcDuqN1tyQbWNnbUH0nZ29zvKiyi5i4PcGtce92pt+ZLa29tY5/U0Tu59LAinxw9X1etGFs1l7FymosI7DMswIr+/7Nj8zuN/celjY/bfbP6Q8LPtCW3phStOWWZAMre2Aje9v63nsVd3Lpoic2huMCK7sBZ7+ueI82Blb8G6fxuy/4nRJZNYFlBJZ51hGsjg9LRY/WdeUX0SN3UUVP9D2XFw+suoXrcu9tblGRA3WNWV90O08JrIruLAqsbFx2uQ/bf+6q+7Y1RmCZRwRW5FxAxYNf4gxL3S4+e1Ip/XC19wVWBNarxynC9sCK3yHV9168bCJ3ctFJuuDzkA76UaE6sBInyIpSI3FGvui+ElhGYBlzTmDVnSKcvyd+FKTjEayKIxkdj2AddKyu+/1TdJyv4rDcsCNVdwmsXqcIW5ZGYBmBZczhjZX4896DSmlgRfby1c8CSVh1zyKKPLNkD+37HKyiI1iR2xB/DlbRbX7FnnsXuTGvqifeHR1YkWekVSda0c2OfOeWfke8jnwOVtHKqisjsMyHNNbmzn3vpMPm+/f2xRUvK4u/M707bn8V4St3mq8Rrbsl7a8ifJWcvqy4zaWfp+VUWnBjiIdXon4St7b9FGHpgpauS/C+at8wigIr/gUaI7DMJxzKurtiTu/1a25FNr+77IVcyd0ILKOxPLwZD73GyhqBZcxJ+1C7TtNlK7IhWVljBJYxxhhjjMAyxhhjjBFYxhhjjDFGYBljjDHGCCxjjDHGGIFljDHGGGMEljHGGGPMhQNrmqavf77q3n7986v6rRo93f09TaVvb7fi33K5XC6Xy612Bdbt3enPQta56VXgcrlXcL8/5t9/C96m2QOD+5l7hFuxTXbZnu/lCiyBxeVyBZb7mSuwBJbAsuPgcgWWwOIKLIElsAQWl8sVWFyBJbAElsDicrkCy/pyBZbAElh2HFyuwBJYXIElsASWwOJyuQKLK7AElsASWFwuV2BZX67AElgCy46DyxVY7meuwLp7YE1/Zv2e9fsFlsDicrkCiyuwBFZ9YC3+ILAEFpfLLQ2s0nE/cwXW5wTWT0LNQ2ovtgSWwEq75b/06fw3O2gulyuwBFb/wNo8UiWwBFZdOd3x622JMw8MXC5XYAksgSWw6oPDDqv7ETsPDFwuV2B9bGDtRZXAek5geeC/45Gzs96sL5crsARWNLAWI7CuH1hHP2TaYXEdseNyBZbA6nYdLK8ivGxgbT5Q2XFwPdfN/czlCqxbBpbrYJ0VWJEf/e04uI6cOSXK5QosV3IXWKlvpM3duh0Hl3uXsBvw9QpKrsASWAIr/o2U2q/ZcXC5z3yxSLp+PMfO9iywBJbAygeWHQeX+5znnG26ridnuxJYAktg9XT/7FnsOLjcB7nBYz9OxcZrzHYlsASWwBJYXC6X2+oWnal0PwssgfW4wPr5OcwOmsvlcju6Tk0KLIElsAQWl8vlunCuwBJYAquf+/7JyY6Sy+VyP/5J/QJLYAksgcXlcrnczm7RmUqBJbAEVqU75no2dpRcLpd7C3fwqyYFlsASWHZYXC6X61WTnV81KbAE1mcG1uLHETssLpfL5R70qkmBJbAElh0Hl8vlcoc+90tgCSyBZcfB5XK53M7PwQq2l8ASWDdwx79axA6Ly+VyuaWnCPfONgosgSWw7LC4XC6X2+05WHUnGQWWwBrhbv5aBjsOLpfL5d70OliR04uPDqxpNpH3CyyBxeVyuVyBFTy9+NDAWkdV+v0Cq87d+62idhxcLpfL/ewruVefXvyoU4QCS2BxuVwuV2AdHToXf/WiwLqHu1dXdhxcLpfL9bsIL/jqxW6B5TlYAovL5XK5Autqp+rOevWiI1g3cBN1ZcfB5XK53PGBVTpXu9Bo9vqoAktg2XFwuVwu9zbuxS802vFXA3kVocCy4+ByuVyuwMq7weRyHax7uOm68g3M5XK5XIHlSu4CS2BxuVwuV2AJLIF1qvvz+lLfwFwul8sVWAJLYAksLpfL5XIFlsC6pPu+PJpvYC6Xy+UKLIElsAQWl8vlcrkCS2Bdz33XlW9gLpfL5QosgXWPwPqdvpTY1tvbrfi3Fe73fxZu45VzS4fL5V7EPfGK29b34u5NA+sh6/vQwGop2QEb9Pzw1TN/9yKXyz30gdD9/NmuI2dXcAWWwLLD4nIFlvtZYAksgSWw7LC4XK7A4godgSWwihZ4UVd2HFwuV2BxhY7AElgCi8vlCiyuwBJYAutKgbWuKzsOLpcrsLhCR2AJLIHF5XIFFldgCSyBdZnA2qwrOw4ulyuwuEJHYAksgcXlcgUWV2AJrB6BNb966Tqk9t4vsOYLvFdXdhxcLldgcYXOEwNrEU/z/937s8ASWFwuV2BxBZbAKjhF+A6p9VErgfUV+NXOdhxcLldgcQWWwMoEllOE6QV+/35nOw4ul7u3v+r+y4DdzwJLYN0msNYhNf/fdWM9ObDeUfXuKj+JcrlcLldgCayaU4QPD6x1VDnUz+VyuVyBJbAEVvEdnY4qgcXlcrlcgSWwml5F+JxThPOo8mRVLpfL5QosgdX/Oljp62N9XGAto0pgcblcLldgCSxXcq8MrL+HrM5ZYDssLpfL5QodgfVpgfU+aiWwuFwulyuwuAKrQ2DNzwkKLC6Xy+UKLK7Aagqs9TOuBBaXy+VyBRZXYNUH1uaT2QUWl8vlcgUWV2BVBtbeSwUFFpfL5XIFFldgFQdW+kIMAovL5XK5AosrsMoCK5FWAovL5XK5AosrsCoC6/dlF9gOi8vlcrlCR2DdLLDeFxEVWFwul8sVWFyB1SGwvtOq5Zc9Cywul8vlCiyBJbC260pgcblcLldgcQVWa2D9nBYs/2XPAovL5XK5AosrsHIHrgQWl8vlcgUWV2C1BtZmXQksLpfL5QosrsASWHZYXC6XK7CEzv0Da5rNZkut33/BwNqrK4HF5XK5XIHFHR1Yi3jabCmBJbC4XC6XK7AEVv0pws3eun5gJepKYHG5XC5XYHEvFFjvPwssgcXlcrlcgSWwigNr8zlYdwmsdF0JLC6Xy+UKLO5VjmBtHsoSWAKLy+VyuQJLYDUF1mKuGVjZuhJYXC6XyxVY3Mu9ivDiR7AEFpfL5XIFlsC633WwrhxYkboSWFwul8sVWFxXchdYdlhcLpcrsISOwDopsIJ1JbC4XC6XK7C4Aktg2WFxuVyuwBI6AuuMwIrXlcDicrlcrsDiCiyBZYfF5XK5AkvoCKwTAqugrgQWl8vlcgUWV2AJLDssLpfLFVhCR2CN3aC/62r9y54FFpfL5XIFlsASWALLDovL5XK5QkdgXSOwvuvqa/XLngUWl8vlcgWWwBJYAssOi8vlcrlCR2BdI7B+6kpg2WFxuVyuwBJYAktg2XFwuVwuV+gIrEsG1ndd/Vr9smeBxeVyuVyBJbAElsCyw+JyuVyu0BFY1wisn7oSWHZYXC6XK7AElsASWHYcXC6XyxU6Tw2saTaR9w8OrHddCSw7LC6XyxVYAusegbWOqvX7N9tLYAksLpfL5QosgRU6Rbg+WLXXYWM26HldCSw7LC6XyxVYAuvzA+v3n98LWPT2/oKDH//9n/ef54FVOm+34t92ce04uFwulzt3z3o84kbcboG191yrvepKN2OvIzqLw1e+gblcLpfL5Q5wRxzB2qwugcXlcrlcLldgVQZW4pjW0YG1risbFpfL5XK53HsEVqKo9upKYHG5XC6XyxVYNdfBWj/na2RgbdaVDYvL5XK5XO4tTxF2uZK7wOJyuVwulyuwrhVYe3Vlw+JyuVwulyuwBBaXy+VyuVyBdYHAStSVDYvL5XK5XK7AKvuCfy7YboG5XC6Xy+UKrD6BlU0rGxaXy+VyuVyBFf2CIweubFhcLpfL5XIFVjSw4mllw+JyuVwulyuwMl9w0YErGxaXy+VyuVyBlf7bmrSyYXG5XC6XyxVYu3VlgblcLpfL5QqsPoH1c+DKAnO5XC6XyxVYfQLrJ60sMJfL5XK5XIHVIbDeB64sMJfL5XK5XIHVIbAWaWWBuVwul8vlCqz6wFofuLLAXC6Xy+VyHxRY02wi7w8E1m8LzOVyuVwu97mBtY6q9Z/LA2uywFwul8vlcp0izETV+n8FFpfL5XK5XIElsLhcLpfL5XIHBtb6uVYCi8vlcrlcrsByBMuGxeVyuVwuV2BZYC6Xy+VyuV5FKLC4XC6Xy+UKrLOugyWwuFwul8vlOkV48C97tsBcLpfL5XIFlsDicrlcLpcrsASWDYvL5XK5XK7AssBcLpfL5XIFlsDicrlcLpcrsAQWl8vlcrlcrsCywFwul8vlcgWWwOJyuVwulyuwBBaXy+VyuVyuwLJhcblcLpfLFVgCi8vlcrlcrsASWFwul8vlcrkCy4bF5XK5XC5XYFlgLpfL5XK5AktgcblcLpfLFVgCi8vlcrlcLvcSgTXNJvJ+gcXlcrlcLldgpQJrHk/zllrHlsDicrlcLpcrsGpOEQosLpfL5XK5Auu0wPrnP197b+m//ee/v+rfEp82+8blcrlcLpcbeOscWEXPwRJYXC6Xy+VyBVYmsLKHrAQWl8vlcrlcgVUQWJvHqAQWl8vlcrlcgVV/mYaKA1oCi8vlcrlcrsDavUzDYjwHi8vlcrlcrsA650ruAovL5XK5XK7AElhcLpfL5XK5AsuGxeVyuVwuV2AJLC6Xy+VyuQJLYHG5XC6Xy+UKLBsWl8vlcrlcgWWBuVwul8vlCqxoYE3lY4G5XC6Xy+UKrJojWBaYy+VyuVyuwBJYXC6Xy+VyuQKLy+VyuVwuV2BZYC6Xy+VyuQJLYHG5XC6XyxVYAovL5XK5XC5XYFlgLpfL5XK5AktgcblcLpfLFVgCi8vlcrlcLnd8YM2vtL4Oqb33W2Aul8vlcrkCazuw5vG0aKm9PwssLpfL5XK5AqvgFOE7pNZHrQQWl8vlcrlcgdUhsJwi5HK5XC6XK7CaAmvvdOG6sQQWl8vlcrlcgZUPrHVCpf/XAnO5XC6XyxVYqcDaPAkosLhcLpfL5Qqs+ss0ZJ+PJbC4XC6Xy+UKrILLNCwmeH0sC8zlcrlcLldguZI7l8vlcrlcrsCyYXG5XC6XyxVYFpjL5XK5XK7AElhcLpfL5XK5AovL5XK5XC5XYFlgLpfL5XK5AktgcblcLpfLFVgCi8vlcrlcLldgWWAul8vlcrkCS2BxuVwul8sVWAKLy+VyuVwuV2BZYC6Xy+VyuQLLAnO5XC6XyxVYAovL5XK5XC5XYNmwuFwul8vlCiwLzOVyuVwu9/GBNc1ms6XW7xdYXC6Xy+VyBdZuYM3jaa+lBBaXy+VyuVyBVX+KcNFSP/8rsLhcLpfL5QqsPoH1/rPA4nK5XC6XK7AqA2vz8JXA4nK5XC6XK7AqA2uvrgQWl8vlcrlcgVUTWHvPbd97gaHA4nK5XC6XK7Ayl2nIXvXKESwul8vlcrkCq+AyDYkjVQKLy+VyuVyuwHIldy6Xy+VyuVyBZYG5XC6Xy+UKLIHF5XK5XC5XYAksLpfL5XK5XIFlgblcLpfL5QosC8zlcrlcLldgCSwul8vlcrlcgWXD4nK5XC6XK7AsMJfL5XK5XIElsLhcLpfL5QosgcXlcrlcLpcrsCwwl8vlcrlcgSWwuFwul8vlCiyBxeVyuVwulyuwLDCXy+VyuVyBJbC4XC6Xy+UKrP3AmmYTeb/A4nK5XC6XK7BSgTWPp3lL7b1fYHG5XC6XyxVYZacI1werNt8vsLhcLpfL5QosgcXlcrlcLpd7RmAF60pgcblcLpfLFVihwIrXlcDicrlcLpcrsPKBVVRXAovL5XK5XK7Ayl+moaiuBBaXy+VyuVyBlblMw2LS7xdYXC6Xy+VyBZYruXO5XC6Xy+UKLAvM5XK5XC5XYFlgLpfL5XK5AktgcblcLpfL5QosGxaXy+VyuVyBZYG5XC6Xy+UKLIHF5XK5XC6XK7C4XC6Xy+VyBZYF5nK5XC6XK7AEFpfL5XK5XIElsLhcLpfL5XIFlgXmcrlcLpcrsAQWl8vlcrlcgSWwuFwul8vlcgWWBeZyuVwulyuwLDCXy+VyuVyBJbC4XC6Xy+VyDwqsaTaR9wssLpfL5XK5AisVWPN4WrTU+88Ci8vlcrlcrsCqP0W4F1Xr/7XAXC6Xy+VyBVbvwPr6qnz79av+rRrlcrlcLpfLjb11DqzN84MCi8vlcrlcrsCqDKxsUQksLpfL5XK5AqsgsDafxi6wuFwul8vlCqz6yzQUPR9LYHG5XC6XyxVYmcs0LCZ6HSwLzOVyuVwuV2B1vpK7BeZyuVwulyuwBBaXy+VyuVyuwLJhcblcLpfLFVgWmMvlcrlcrsASWFwul8vlcrkCi8vlcrlcLldgWWAul8vlcrkCS2BxuVwul8sVWAKLy+VyuVwuV2BZYC6Xy+VyuQJLYHG5XC6XyxVYAovL5XK5XC5XYFlgLpfL5XK5AssCc7lcLpfLFVgCi8vlcrlcLldg2bC4XC6Xy+XeK7Cmv7NZUXt/K7C4XC6Xy+UKrMwRrM2ESv+vBeZyuVwulyuwBBaXy+VyuVyuwLJhcblcLpfL/aTA8hwsLpfL5XK5AssRLC6Xy+VyuVyBZYG5XC6Xy+UKLIHF5XK5XC5XYBVcB2v9XCvPweJyuVwulyuwXMmdy+VyuVwuV2BZYC6Xy+VyuQJLYHG5XC6XyxVYAovL5XK5XC5XYFlgLpfL5XK5AssCc7lcLpfLFVgCi8vlcrlcLldg2bC4XC6Xy+UKLAvM5XK5XC5XYAksLpfL5XK5XIHF5XK5XC6XK7AsMJfL5XK5XIElsLhcLpfL5QosgcXlcrlcLpcrsCwwl8vlcrlcgSWwuFwul8vlCiyBxeVyuVwulzsssKa/sxdSm38rsLhcLpfL5QqszBGszcCav3PxAQKLy+VyuVyuwCoOrL1jWgKLy+VyuVyuwGoKLKcIuVwul8vlCqyegfV+57qxBBaXy+VyuVyB1XqKUGBxuVwul8sVWAKLy+VyuVwu98KvInSKkMvlcrlcrsCquQ7WZkV5kjuXy+VyuVyB5UruXC6Xy+VyuQLLAnO5XC6XyxVYAovL5XK5XK7AElhcLpfL5XK5AssCc7lcLpfLFVgWmMvlcrlcrsASWFwul8vlcrkCy4bF5XK5XC5XYFlgLpfL5XK5AktgcblcLpfLFVgCi8vlcrlcLldgWWAul8vlcrkCS2BxuVwul8sVWAKLy+VyuVwuV2BZYC6Xy+VyuQJLYHG5XC6XyxVYO4E1/Z1ES63/VmBxuVwul8sVWJkjWAKLy+VyuVwud1Bg/bxfYHG5XC6XyxVYfQLr/U6BxeVyuVwuV2AJLC6Xy+VyudzrBdb8PQKLy+VyuVyuwOoTWIsRWFwul8vlcgVWn1cROoLF5XK5XC5XYFVeBytxNSyBxeVyuVwuV2C5kjuXy+VyuVyuwLLAXC6Xy+VyBZYF5nK5XC6XK7AEFpfL5XK5XK7AsmFxuVwul8sVWBaYy+VyuVyuwBJYXC6Xy+VyuQKLy+VyuVwuV2BZYC6Xy+VyuQJLYHG5XC6XyxVYAovL5XK5XC5XYFlgLpfL5XK5AktgcblcLpfLFVgCi8vlcrlcLldgWWAul8vlcrkCywJzuVwul8sVWKvAmv7OZkXt/a3A4nK5XC6XK7AyR7A2E2pRWgKLy+VyuVyuwGoKrPQHCCwul8vlcrkCS2BxuVwul8vlXimwPAeLy+VyuVyuwOoZWHvPf7fAXC6Xy+VyBVZNYCXeb4G5XC6Xy+UKrOLASh/WssBcLpfL5XIFVuY6WItrMUyrEVhcLpfL5XIFliu5c7lcLpfL5QosC8zlcrlcLldgWWAul8vlcrkCS2BxuVwul8vlCiwbFpfL5XK5XIFlgblcLpfL5QosgcXlcrlcLpcrsLhcLpfL5XIFlgXmcrlcLpcrsAQWl8vlcrlcgSWwuFwul8vlcgWWBeZyuVwulyuwBBaXy+VyuVyBJbC4XC6Xy+VyBZYF5nK5XC6XK7AsMJfL5XK5XIElsLhcLpfL5XKPC6zp72xW1N7fCiwul8vlcrkCK3MEay+w9v5WYHG5XC6XyxVYxYG1eM/6fy0wl8vlcrlcgSWwuFwul8vlcgWWDYvL5XK5XK7AssBcLpfL5XIFlsDicrlcLpfL9SpCLpfL5XK53LtcB2t9vSvXweJyuVwulyuwXMmdy+VyuVwuV2BZYC6Xy+VyuQJLYHG5XC6XyxVYAovL5XK5XC5XYFlgLpfL5XK5AssCc7lcLpfLFVgCi8vlcrlcLldg2bC4XC6Xy+UKLAvM5XK5XC5XYAksLpfL5XK5XIHF5XK5XC6XK7AsMJfL5XK5XIElsLhcLpfL5QosgcXlcrlcLpcrsCwwl8vlcrlcgSWwuFwul8vlCqxkYE2zEVhcLpfL5XIFVmtgLaIq2FgCi8vlcrlcrsASWFwul8vlcrkCywJzuVwul8v1HCwLzOVyuVwuV2A5gsXlcrlcLpcrsGxYXC6Xy+VyBZYF5nK5XC6X6zlYnoPF5XK5XC5XYLmSO5fL5XK5XK7AssBcLpfL5XIFlsDicrlcLpcrsAQWl8vlcrlcrsCywFwul8vlcgWWBeZyuVwulyuwBBaXy+VyuVyuwLJhcblcLpfLFVgWmMvlcrlcrsASWFwul8vlcrkCi8vlcrlcLldgWWAul8vlcrkCS2BxuVwul8sVWAKLy+VyuVwuV2BZYC6Xy+VyuQJLYHG5XC6XyxVY8cCa/o7A4nK5XC6XK7A6BNa8q4KNJbC4XC6Xy+UKrN3Aih+1ElhcLpfL5XIFVkFgOUXI5XK5XC6X2zOw3l0VbyyBxeVyuVwuV2BFTxEKLC6Xy+VyuQJLYHG5XC6Xy+Ve+FWEThFyuVwul8vldgssT3LncrlcLpfLdSV3LpfL5XK5XIFlgblcLpfL5QosgcXlcrlcLldgCSwul8vlcrlcgWWBuVwul8vlCiyBxeVyuVwuV2AJLC6Xy+VyuVyBZYG5XC6Xy+UKLAvM5XK5XC5XYAksLpfL5XK5XIFlw+JyuVwulyuwLDCXy+VyuVyBJbC4XC6Xy+VyBRaXy+VyuVyuwLLAXC6Xy+VyBZbA4nK5XC6XK7AEFpfL5XK5XO5FAmv6MwKLy+VyuVyuwBJYXC6Xy+VyuZcMrJ+0ElhcLpfL5XK5fQLr3VUCi8vlcrlcLldgcblcLpfL5V4vsOZRJbC4XC6Xy+Vy+wTWYgQWl8vlcrlcgdXtOliOYHG5XC6Xy+UKLC6Xy+VyuVxXcrfAXC6Xy+VyBZbA4nK5XC6XK7AEFpfL5XK5XK7AssBcLpfL5XIFlsDicrlcLpcrsAQWl8vlcrlcrsCywFwul8vlcgWWBeZyuVwulyuwBBaXy+VyuVyuwLJhcblcLpfLFVgWmMvlcrlcrsASWFwul8vlcgWWwOJyuVwul8sVWBaYy+VyuVyuwBJYXC6Xy+VyBZbA4nK5XC6Xyz0xsKbZCCwul8vlcrkCqzWw5lEVbyyBxeVyuVwuV2BFTxEKLC6Xy+VyuQJLYHG5XC6Xy+VeOLA8B4vL5XK5XC63Z2DF60pgcblcLpfLFVj5wCqqK4HF5XK5XC5XYOUv0+A6WFwul8vlcrk9L9OwGIHF5XK5XC5XYLmSO5fL5XK5XK7AssBcLpfL5XIFlsDicrlcLpcrsAQWl8vlcrlcrsCywFwul8vlcgWWBeZyuVwulyuwBBaXy+VyuVyuwLJhcblcLpfLFVgWmMvlcrlcrsASWFwul8vlcrkCi8vlcrlcLldgWWAul8vlcrkCS2BxuVwul8sVWAKLy+VyuVwuV2BZYC6Xy+VyuQJLYHG5XC6XyxVYscCaZiOwuFwul8vlCqw+gbX4g8DicrlcLpcrsOoDaxFVwcYSWFwul8vlcgVW/8AyxhhjjHngHBhYxhhjjDFGYBljjDHGCCxjjDHGmJsG1qvqVYTGGGOMMSYTWKXXwTLGGGOMManA6vZJ2yqtpfPaG7Hi38ZfU9D3Hit6OcPV7ufs1x750hJ/e/RPC3ufP+KmP+bon3AWn79oE6r7eo/YQqKv4hl+P+99VxZtz+sPe+z2XPo9ftB+o+j+P+LzV+/rDrrNV7ifD32U6XJLjt2PN/6r0q+ty9PI6kql1x3VuCcac18tHjYO2lqKrnY74BYG772sm73nxwdWY5wddz8nvjsibvrLHHOovuX+WcTZodvz3uePbM/pj+n1Y1jwvm3ZbyQ+Z91VuOOPbsHPX/Q5K+7/7H64+/0weDuPP8p0uSWXC6yOn6R6jzYysI4LlEMDa8BCH/ENPPhBNPiDY/sW2LKFl97PiX/b/WZXR2rwgWrYhtH4HTr4fi7antMf0+t+Pv2Bv+L+jxwPLv388c9Zd/9X3+bB98PRjzK9bonA6lCsjefpXs0nVT81sOKngT4jsIYdDao4RRg5WnDQDnHz9ETwgX/v4P/gDaNiL/EBgdXxfi49CtJ3v3HTwDo6hqrv50sF1t4J0M8PrJHPZ+q1O6g4rbl5gHTkHVUXlI23ufR4Q9FO54hbGNnSinaU4x/4i073xOPsuPt57/MHH/j3btuYwMquflHxPPx+rvshqmK/kYiJ0vs/vsS9AuuIo3qRR9XS+zmx/zx6O49/F3x4YA3us15PhDrlaFBLTTYetHuNOlS5/skpsm8dfCYo+1N13ZGhI46mFN3P6QeD4+7n0iMKiRsz7H4ufSxMp/aw7Xnz8we35+y/PSuwSvcb6aNxpfd/vDni+426z3nEEazG/XPH+/mgR5lPDqzxpxcbX5F3YmB1abLrL1bH5yWMD6yiB9qDXiwW/PzZl/yU/tuRgVX0oHX0i/Jszx3v5/bAGvM56/55+xGs7Mmv7qc1B2zPAuuQvUCvfzjsBrS88vG1cyD6yoF1xDGAilfHxP/t0dt231cRvs578nXwSg3D7ufEd0fRqwizhw8H7P3usj1XvIpwwCNl3dHfilNXkc/Z5TYHH9pbPmeX2xw8EtzrVOzI/farx+s6TwisXhdnGnxtp8bDSKdcU6rLBWZabvOhW0vpdYMOuoVF9177dYMGB9ar/Bo/7ufqn8jjLyY47juu6PMX3c9H76aC20CX/UbLdh65mFnpdc4qPmfpocfb3c/nbkUnH8EyxhhjjHnmCCxjjDHGGIFljDHGGCOwjDHGGGMEljHGGGOMEVjGGGOMMQLLGGOMMUZgGWOMMcYYgWWMMcYYI7CMMcYYYwSWMcZ80A4x+Us/0h+T+JxH/+4jY4zAMsac3xDtv2ar+6+d7/5JugdWy80TWMYILGPM5wdWx8f+D0sHgWWMEVjGmM6BFf8V9z9/Xn989rMt/mH6FsY/eHEb0l9F4v2bSvDj059HYxkjsIwxTwysdaNsNtPeH+KfLXs0qPSD15WT/ioitzMbWEVfpsAyRmAZYz4/sBLHlhItFfnb6ve3f3D8nxxE1N0SY4zAMsZ8SGBF3r958usWgVV0Ck9gGWMEljHmnMBKxMTFj2A1VpHAMsYILGNMa2CNeQ7W0YGVPdIWf+LUplLxHCx1ZYzAMsY8NLBe5a8iXP9V8FWExwVW8KuIv/ov/XUFP4/AMkZgGWOM2U7PunJ1JXdjBJYxxpj6CHMnGGMEljHGCCxjjMAyxhhjjBFYxhhjjDHPnP8Bo93IuxAbHUoAAAAASUVORK5C)

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

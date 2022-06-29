# syntax=docker/dockerfile:1
FROM r-base:4.1.0
MAINTAINER xwj
RUN R -e "install.packages('devtools', repos=c('https://mirrors.tuna.tsinghua.edu.cn/CRAN/', 'http://cran.rstudio.com/'))" 
CMD ['date']
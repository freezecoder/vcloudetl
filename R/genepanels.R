#Gene panels code

library(tidyverse)

getPanel<-function(file="presets.yaml"){
  return(suppressWarnings(yaml.load_file(file)))
}

panelNames<-function(s){
  names(s$panels)
}

panelElements<-function(s,verb="genes",id="Lung1"){
  
  if (verb %in% names(s$panels[[id]])){
    s$panels %>% map(`[`,verb) %>% as_tibble() %>% select(id) %>% unnest() %>% unlist() %>% unname()
    
  }else{
    return(NULL)
  }
  
}

panelAllElements<-function(s,verb="genes",panelname="Lung1"){
  panelNames(s) %>% map(~panelElements(s,id=.x,verb=verb)) %>% unlist() %>% unique()
}

#Apply  panel filter
filterByPanel<-function(indata="",paneldata=getPanel(),panelname="Lung1",testing=F){
  printf("Filtering using %s panel",panelname)
  ldat=NULL
  tdat=NULL
  gdat=NULL
  dat<-copy(indata)
    setnames(dat,tolower(names(dat)))
  setkey(dat,"hugo_gene_symbol")
  #filter by genes
  igenes=panelElements(paneldata,verb="genes",id=panelname)
  itrans=panelElements(paneldata,verb="transcripts",id=panelname)
  ilocation=panelElements(paneldata,verb="location",id=panelname)

  dat$start_position=as.double(dat$start_position)

  #Make the HGVSC column  
  if ("hgvsc" %in% names(dat)){
    dat$transcriptname=gsub(":","",str_extract(dat$hgvsc,"(\\S)+:"))
    dat$cnomen=gsub(":","",str_extract(dat$hgvsc,":(\\S)+"))
    
    #Transcript matching
    if (!is.null(itrans)){
      itrans=gsub("\\.\\d+$","",itrans)
     # print(itrans)
      pt=paste(itrans,collapse="|")
      tdat=subset(dat,grepl(pt,hgvsc,perl=T))
      #print(dat$transcriptname)
      printf("%s Transcripts matching",nrow(tdat))
    }
  
  }
  #Gene symbol matching
   if (!is.null(igenes)){
       gdat=subset(dat,hugo_gene_symbol %in% igenes)
       printf("%s Gene symbols matching",nrow(gdat))
   }
  if (!is.null(ilocation)) {
    ilocation %>% map(~ coordMatch(dat=dat,loc=.x)) %>% tibble() %>% unnest() -> ldat
  ldat=as.data.table(ldat)
    printf("Range Location Matches: %s",nrow(ldat))
  #  print("Not running location filter yet")
  }

  #Collect results together from transcript,gene and location matches
  dat=rbind(tdat,gdat)
  if (!is.null(ldat)){
    dat=rbind(dat,ldat)
  }
  dat=unique(dat)
  
  if (testing==TRUE){
  #Testing
  dat=subset(dat,select=c("hugo_gene_symbol","symbol","location","allele","transcriptname","hgvsp"))
  dat=unique(dat)
  }
  dat=as.data.table(dat)
  dat=dat[hugo_gene_symbol !="." | hugo_gene_symbol !="-" | hugo_gene_symbol !=""]
  return(dat)
  
}


coordMatch<-function(dat="",loc=""){
  coord=unlist(str_split(loc,pattern=":"))
  ps=unlist(str_split(coord[2],pattern="\\-"))
  pstart=as.numeric(ps[1])
  pend=as.numeric(ps[2])
  pchr=coord[1]
  mres=subset(dat, as.numeric(start_position) >= pstart & as.numeric(start_position) <= pend & chrom==pchr )
  mres=as.data.frame(mres)
  #if (nrow(mres)>0){
  return(mres)

  
}
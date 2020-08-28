## Variant ETL -> from TSV format
## Loading  with elastic from input files residing on the cloud
## add unique keys to each record before indexing in ES
## Get methods

library(elasticsearchr)
#setCloudEnv()

printf<-function(...){
  print(sprintf(...))
}
#' Title
#'
#' @return
#' @export
#'
#' @examples
deleteVariants<-function(){
  elastic("http://localhost:9200", "variants", "vdata") %delete% TRUE
}

testElasticLoad<-function(){
  library(elasticsearchr)
  es <- elastic("http://localhost:9200", "variants", "vdata")
  
  
}


# Title
#
# @export
# @example
hashedResultsdb<-function(filter=NULL,uid=NULL,pattern=NULL){
	#get s3 list of files in a data frame
	return(df)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
getSourceVariantFiles<-function(filepattern="variantTable.tsv$",uid=19304){
  mydf=hashedResultsdb(filter=filepattern,uid=uid)
  mydf$prefixpath=paste0(mydf$subpath,"/",mydf$filebasename)
  mydf=unique(mydf,by="sampleid")
  return(mydf)
  
}
################################################################################
################################################################################
################################################################################
################################################################################
#Main ETL Process----

#Load multiple variant files ======
loadVariantFiles<-function(mydf=getSourceVariantFiles()){
  print(names(mydf))
  esrc <- elastic("http://localhost:9200", "loadstat", "data")
  #Check what's been loaded
  for_everything <- query('{
  "match_all": {}
  }')
  esrc %search% for_everything %>% filter(type=="SNV") %>% select(projectid,sampleid) %>% distinct -> loaded
  loaded$isloaded=1
  print(nrow(mydf))
  #left_join(as_tibble(mydf),loaded) -> mydf
  #mydf = subset(mydf, sampleid %in% loaded$sampleid)
  mydf=merge(mydf,loaded,by=c("projectid","sampleid"),all.x=T)
  print(table(mydf$isloaded))
  mydf=subset(mydf, is.na(isloaded))
  #mydf=as.data.frame(mydf)
  
  if (nrow(mydf)==0){
    printf("%s Files already in database store",nrow(loaded))
    print("No files Loaded")
    return("No Files Loaded")
  }else{
    printf("%s Files to load..",nrow(mydf))
  }
    for (idx in 1:nrow(mydf)){
       #print(idx)
       sfile=mydf[idx]$prefixpath
       #print(names(mydf[idx]))
       if (is.null(file))
         next
       bucket=mydf[idx]$Bucket
       #print(bucket)
       bs=mydf[idx]$filebasename
       #printf("Loading %s %s\n",bucket,bs)
       sampleid=mydf[idx]$sampleid
       sampleid=gsub("\\.|\\s+","",sampleid)
       sampleid=tolower(sampleid)
       pjid=mydf[idx]$projectid
       #printf("Load %s %s %s",file,bucket,pjid)
       printf("Load %s %s %s",sfile,bucket,pjid)
       #print(sfile)
       tryCatch({
         
        rep=loadVariantBulk(s3path=sfile,s3bucket = bucket,pid=pjid,doFilter = FALSE)
       },
       error=function(e){
        print(e) 
       }
       )
       print(rep)  
       loadtime=Sys.Date()
       created=mydf[idx]$time
       
       
       dfstatus=data.frame(
         type="SNV",
         id=paste0(as.numeric(as.POSIXct("2013-09-16 2:13:46 EST")),sample(234243,1)),
         created=created,
         loaded=loadtime,
         sampleid=mydf[idx]$sampleid,
         projectid=pjid,
         nrows=rep$inrecords,
         loaded=rep$loaded
       )

      esrc %index% dfstatus
    }
  return("Files have been loaded")
  
  
}

#main method for loading a file from an s3 bucket location or from a local file 
#S3 inputs require access credentials via the aws.s3 package
#' Insert data from a file/cloud path into ES store
#'
#' @param mydf 
#' @param pid 
#' @param tblclass 
#' @param istart 
#' @param tblname 
#' @param infile 
#' @param source 
#' @param s3path 
#' @param s3bucket 
#' @param doFilter 
#'
#' @return
#' @export
#'
#' @examples
loadVariantBulk<-function(pid="pjx",tblclass="vdata",skipLoad=FALSE,istart=1,tblname="snvs",infile=NULL,source="cloud",
                          s3path=mydf$prefixpath[1],s3bucket=NULL,doFilter=TRUE){
  es <- elastic("http://localhost:9200", tblname, tblclass)
  #get data
  print(s3bucket)
  if (source=="cloud"){
    print(s3path)
    dres=s3read_using(FUN=fread,object=s3path,bucket=s3bucket)
  }else {
    dres=fread(infile)
  }
  #Remove duplicate column names
  dres$TYPE=NULL
  dres$CHROM=NULL
  dres$entrez_gene_id=NULL
  
  nms=tolower(gsub("\\s+","_",names(dres)))
  nms=gsub("#","",nms)
  nms=cleaned_field_names(nms)
  
  setnames(dres,nms)
  #Filter out all synonymous variants ####
  if (doFilter==TRUE){
    dres=subset(dres,consequence !="synonymous_variant")
    dres=subset(dres,canonical=="YES")
  }
  
  if (nrow(dres)==0){
    print("No variants to load")
    return(NULL)
    
  }
  dres=unique(dres)
  inrecords=nrow(dres)
  dres$projectid=pid
  dres$rcnt=1:nrow(dres)
  #important, append a unique ID as a new field in the input table
  dres[["id"]]=gsub("\\-|\\_|\\/","",paste0(dres$rcnt,dres$sample_id,dres$uploaded_variation))
  dres$rcnt=NULL
  
  if (skipLoad != FALSE){
    index_bulk_dataframe(es,dres)
    printf("Indexed %s records",nrow(dres))
  }
  #Load the summary tables
  loadVariantSummaries(dres)

  return(list(
    inrecords=inrecords,
    loaded=nrow(dres)
    
  ))
    
}

#Load summaries of the variants from input data frame
loadVariantSummaries<-function(dres,tblclass="vdata",tblname="snvs"){
  #Add some aggregated data
  tryCatch({
      genesummary= dres[hugo_gene_symbol!="-",by=c("sample_name","hugo_gene_symbol"),list(records=length(type), clinsig=paste(unique(clin_sig),collapse=","))]  
      genesummary[["id"]]=paste0(genesummary$sample_name,sample(1:3030,1),1:nrow(genesummary))
      tblname="genesum"
      es2 <- elastic("http://localhost:9200", tblname, tblclass)
      index_bulk_dataframe(es2,genesummary)
  },error=function(e){
    print("Error with gene summary")
  })
  tryCatch({
  samplegene=dcast.data.table(dres[hugo_gene_symbol!="-"],sample_name ~hugo_gene_symbol)
  samplegene[["id"]]=paste0(samplegene$sample_name,sample(1:3344,1),1:nrow(samplegene))
  tblname="smgene"
  es3 <- elastic("http://localhost:9200", tblname, tblclass)
  index_bulk_dataframe(es3,samplegene)
  },error=function(e){
    print("Error with Sample  summary calculation")
  })
  
}

loadBigIndels<-function(pid="pjx",smid="abc123",annofile="/udata/workspace/refdata/gencode.v19.annotation.genes.gtf",tblclass="vdata",skipLoad=FALSE,istart=1,tblname="bigdel",infile=NULL,source="cloud",s3path=mydf$prefixpath[1],s3bucket=mydf$Bucket[1],doFilter=F){
  es <- elastic("http://localhost:9200", tblname, tblclass)
  #get data
  print(tblname)
  if (source=="cloud"){
    print(s3path)
    dres=s3read_using(FUN=fread,object=s3path,bucket=s3bucket)
  }else {
    dres=fread(infile)
  }
  
  setnames(dres,c("chrom","start","end","name","reads","readids"))
  dres[,size:=abs(start-end)]
  
  dres$sample_id=smid
  dres$projectid=pid
  
  if (nrow(dres)==0){
    print("No variants to load")
    return(NULL)
    
  }
  dres=unique(dres)
  inrecords=nrow(dres)
  dres$projectid=pid
  dres$rcnt=1:nrow(dres)
  #important, append a unique ID as a new field in the input table
  dres[["id"]]=gsub("\\-|\\_|\\/","",paste0(dres$rcnt,dres$name))
  dres$rcnt=NULL
  
  index_bulk_dataframe(es,dres)
  printf("Indexed %s records",nrow(dres))

  #Index annotated indels ####
  if (file.exists(annofile)){
    annotbl=paste0(tblname,"anno")
    print("Running intersection on known genes")
    printf("Loading into %s",annotbl)
    dres$readids=NULL
    annodf=mapIndelsToAnnotations(df=dres,annofile=annofile)
    if (!is.null(annodf)){
      printf("Annotated =  %s",nrow(annodf))
      annodf=as.data.table(annodf)
      setnames(annodf,"id","deletionid")
      annodf$rcnt=1:nrow(annodf)
      annodf[["id"]]=gsub(":|\\-|\\_|\\/","",paste0(annodf$rcnt,annodf$deletionid))
      annodf$rcnt=NULL
      esa=elastic("http://localhost:9200", annotbl, tblclass)
      index_bulk_dataframe(esa,annodf)
    }
  }
  #Return final count ####
  
  return(list(
    inrecords=inrecords,
    loaded=nrow(dres)
    
  ))
  
}

#Intersect to map indels to genes----
# Run before loading into ES
mapIndelsToAnnotations<-function(df,annofile="/udata/workspace/refdata/gencode.v19.annotation.genes.gtf",minsupport=1){
  len=length(names(df))
  df=subset(df,reads>minsupport)
  if (nrow(df)==0){
    print("warning:empty input")
    return(NULL)
  }
  tmpbed=tempfile(pattern = "bigindel.", tmpdir = tempdir(), fileext = ".bed")
  write.table(df,tmpbed,col.names = F,quote=F,sep="\t",row.names = F)
  cmd=sprintf("intersectBed -a %s -b %s -wa -wb",tmpbed,annofile)
  res=fread(cmd=cmd)
  if (nrow(res)>0){
    nms=names(res)
    nms[1:len]=names(df)
    setnames(res,nms)
    #uniquify results
    res=unique(as.data.table(res),by=c("chrom","start","end","name","reads","projectid","sample_id"))
    return(subset(res,select=c(1:len,16:length(names(res)))))
  }else{
    print("No result")
    return(NULL)
  }
}


loadCloudIndels<-function(force=FALSE){
  print("Finding cloud files matching large indels...")
  mydf=hashedResultsdb(filter="deletions.bed$",uid=19304)
  print(nrow(mydf))
  for_everything <- query('{
  "match_all": {}
  }')
  mydf$prefixpath=paste0(mydf$subpath,"/",mydf$filebasename)
  
  if (force==FALSE){
    esrc <- elastic("http://localhost:9200", "loadstat", "data")
    esrc %search% for_everything %>% filter(type=="BIGDEL") %>% select(projectid,sampleid) %>% distinct -> loaded
    loaded$isloaded=1
    mydf=merge(mydf,loaded,by=c("projectid","sampleid"),all.x=T)
    print(table(mydf$isloaded))
    mydf=subset(mydf, is.na(isloaded))
  }
  
  if (nrow(mydf)==0){
    printf("%s Big deletion files loaded",nrow(loaded))
    print("No new files to load at this time")
    return("No new files to load")
  }
  
  lapply(1:nrow(mydf),function(x){
    print(mydf[x]$filebasename)
    sid=mydf[x]$sampleid
    pjid=mydf[x]$projectid
   tryCatch({
     rep=loadBigIndels(pid=pjid,smid=sid,s3path = mydf[x]$prefixpath,s3bucket = mydf[x]$Bucket)
    
    loadtime=Sys.Date()
    created=mydf[x]$time
    
    dfstatus=data.frame(
      type="BIGDEL",
      id=paste0(as.numeric(as.POSIXct(Sys.time())),sample(234243,1)),
      created=created,
      loaded=loadtime,
      sampleid=sid,
      projectid=pjid,
      nrows=rep$inrecords,
      loaded=rep$loaded
    )
    
    esrc %index% dfstatus
   },error=function(e){
     printf("Failed to load %s: %s",x,sid)
     print(e)
   })
    
  })
  
}

loadMetaData<-function(){
  
  
}

################################################################################
################################################################################
################################################################################
################################################################################
#Query variant store in database -----

# This is a base testing method to return a dataframe
# Main query method
#' Retrieve variants from the database
#'
#' @param tblname 
#' @param tblclass 
#' @param qpjid 
#' @param qgene 
#' @param qsample 
#' @param qtype 
#' @param all 
#'
#' @return
#' @export
#'
#' @examples
getVariants<-function(tblname="snvs",tblclass="vdata",everything=F,qpjid="0C734366DE0C36BBE56C1CF2FB74C24A",qgene="TP53",qsample="048_EK_8-308_S42",qtype="sample",allcols=TRUE){
 
  selected_fields <- select_fields('{
  "includes": ["hugo_gene_symbol","feature","location","chromosome","ref","alt","consequence","type","sample_id","canonical","hgvsp","impact","existing_variation","dp","allele_fraction"]
   }')
  
  for_everything <- query('{
  "match_all": {}
  }')
  es <- elastic("http://localhost:9200", tblname, tblclass)
  
  if (allcols ==TRUE){
    selected_fields <- select_fields('{
  "excludes": ["*_dp","*_DP"]
    }')
  }
  
  #Match by gene
 #gq=list(match=list(hugo_gene_symbol=qgene))
  gq=list(term=list(hugo_gene_symbol.keyword=list(value=qgene)))
  match_gene<-query(toJSON(gq,auto_unbox = T))
  

  #Match by project ID

  #Return everything
  if (everything==TRUE){
    es  %search% (for_everything + selected_fields) %>% unique() ->df
    printf("All data returned %s rows",nrow(df))
  }else {
  
    #Gene query
    if (qtype == "gene"){
      printf("Gene Query= %s",qgene)
      es  %search% (match_gene + selected_fields) ->df
      print(nrow(df))
    }else if (qtype=="sample"){
      sq=list(term=list(sample_id.keyword=list(value=qsample)))
      
      match_sample<-query(toJSON(sq,auto_unbox = T))
      print(toJSON(sq,auto_unbox = T))
      printf("Sample Query=%s",qsample)
      #Sample query
      es  %search% (match_sample + selected_fields) ->df
     print(nrow(df))
    }else if (qtype=="project"){
      sq=list(match=list(projectid=qpjid))
      match_project<-query(toJSON(sq,auto_unbox = T))
      printf("Querying project = %s",qpjid)
      es  %search% (match_project + selected_fields) ->df
      
    }
  }
  
  return(df)
  
}

#Returns sample IDs in the noSQL store ====
getSampleIds<-function(tblname="snvs",tblclass="vdata"){
  
  selected_fields <- select_fields('{
  "includes": ["sample_id"]
   }')
  es <- elastic("http://localhost:9200", tblname, tblclass)
  
  es  %search% (match_sample + selected_fields) 
  uquery= aggs('{
    "langs" : {
        "terms" : { "field" : "sample_id.keyword",  
          "size" : 99999 }
      }
    }')
  es  %search% uquery ->df
  names(df)=c("sample_id","doc_count")
  return(df)
  
}

#Return DF of samples as rows and genes as columns ====
getSampleGeneCounts<-function(tblname="smgene",tblclass="vdata"){
  selected_fields <- select_fields('{
  "excludes": ["id"]
   }')
  by_smname <- sort_on('{"sample_name.keyword": {"order": "asc"}}')
  
  elastic("http://localhost:9200", tblname, tblclass) %search% (for_everything + selected_fields + by_smname)  -> smgene
  
  cols=grep("sample_name",names(smgene),invert=T,value=T)
  smgene=subset(smgene,select=c("sample_name",cols))
  smgene=as.data.table(smgene)
  return(smgene)
}

#Return all files that were loaded
getLoadedVariantFileStatus<-function(){
  esrc <- elastic("http://localhost:9200", "loadstat", "data")

    
    for_everything <- query('{
  "match_all": {}
  }')
    
    esrc %search%  for_everything -> rdf
    return(rdf)
  
}

#Return clinsig summary,example of a summed up table
getClinSigSummary<-function(){
  elastic("http://localhost:9200", "genesum", "vdata") %search% for_everything ->tt
  tt=as.data.table(tt)
  tt=tt[clinsig!="-"]
  fun=function(x){return(paste(unique(x),collapse=","))}
  newt=dcast.data.table(tt,sample_name ~ hugo_gene_symbol,fun.aggregate = fun, value.var = "clinsig") 
  return(newt)
}


#Return the variant object from a db Call
#' Make a SNV variant object from the returned variant data frame
#'
#' @param variantdf 
#'
#' @return
#' @export
#'
#' @examples
getSNVobject<-function(variantdf){
  #Get the snvObject code
  source("protocols/mergevep/snvObject.R")
  sampledat=subset(variantdf,select=sample_id)
  sampledat=unique(sampledat)
  td=snvObject$new(variantData = variantdf,featureData = NULL,sampleData = sampledat)
  return(td)
  
}

getDeletions<-function(sid=NULL,pid=NULL,index="bigdel",class="vdata",brief=T){
  for_everything <- query('{
  "match_all": {}
  }')
  elastic("http://localhost:9200", index, class) %search% for_everything ->res
  
  if (!is.null(sid)){
    res=subset(res,sample_id %in% sid)
  }
  if (!is.null(pid)){
    res=subset(res,projectid %in% pid)
  }
  if (brief==TRUE)
     res$readids=NULL
  
  res=unique(res)  
  
  return(res)
}

#Return annotatd deletions
getAnnoDeletions<-function(sid=NULL,index="bigdelanno",class="vdata"){
  es<-elastic("http://localhost:9200", index,class)
  sq=list(term=list(sample_id.keyword=list(value=sid)))
  
  match_sample<-query(toJSON(sq,auto_unbox = T))
  printf("Annodel: Sample Query=%s",sid)
  #Sample query
  es  %search% (match_sample) ->df
  itot=nrow(df)
  
  df=as.data.table(df)
  df=unique(as.data.table(df),by=c(1:9))
  printf("found %s records (%s)",nrow(df),itot)
  
  df$genesymbol=stringr::str_extract(df$V18,"gene_name\\s+(\\S+)")
  setnames(df,"V18","GencodeV19Annotation")
  df$genesymbol=gsub(";|gene_name|\\\"","",df$genesymbol)

  return(df)
  
}

deleteDeletionsIndex<-function(){
  elastic("http://localhost:9200", "bigdel", "vdata") %delete% TRUE
  elastic("http://localhost:9200", "bigdelanno", "vdata") %delete% TRUE
}

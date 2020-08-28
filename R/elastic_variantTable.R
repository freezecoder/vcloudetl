
#Generic table widget for looking at variants in detail
library(shiny)
library(yaml)
library(data.table)
library(stringr)
library(shinyWidgets)
library(DT)
library(shinyBS)
library(countup)
library(shinyjqui)

source("genepanels.R")
#ES functions for loading data
source("variantETL.R")

### General function ----------
shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))
  }
  inputs
}

printf<-function(...){
  print(sprintf(...))
}

xbutton<-function(){
  tags$button(
    type = "button",
    class = "close",
    'data-dismiss' = "modal",
    'aria-hidden' = "true",
    icon("times")
  )
}

joinCustomDb<-function(dat="",xldat="",showall=T){
  #X-col names
  xkey=c("chromosome","hugo_gene_symbol","cnomen","transcriptname")
  #Y-column names
  ykey=c("jChrom","jGene","jcNomen","jTranscript")
  #Make the join columns as copy as they will disappear when join happens
  xldat$jChrom=xldat$Chrom
  #xldat$jMUT_START=xldat$MUT_START
  xldat$jGene=xldat$Gene
  xldat$jcNomen=xldat$cNomen
  xldat$jTranscript=xldat$Transcript
  #Try HGVSC 
  if ("hgvsc" %in% names(dat) && nrow(xldat)>0){
    showNotification("Join is possible",type="warning",duration=2)
    dat$transcriptname=gsub(":","",str_extract(dat$hgvsc,"(\\S)+:"))
    dat$start_position=as.double(dat$start_position)
    dat$cnomen=gsub(":","",str_extract(dat$hgvsc,":(\\S)+"))
    #Set columns in data used for the join
    dat=merge(dat,xldat,by.x=xkey,by.y=ykey,all.x=T)
    #dat=dat[order(cNomen)]
    return(dat) 
  }else{
    print("No Join Happned")
    return(dat)
  }
  
}

mapToClinvar<-function(dat="",file="clinvar.rds"){
  
  if (file.exists("clinvar.rds")){
      clinvar=readRDS(file)
      clinvar=as.data.table(clinvar)
      clinvar=unique(clinvar,by="VariationID")
      
      clinvar$Chromosome=paste0("chr",clinvar$Chromosome)  
    
      
      if (nrow(dat) >0){
        ycols=c("Chromosome","Start","ReferenceAllele","AlternateAllele")
        xcols=c("chromosome","start_position","ref","alt")
        dat$start_position=as.numeric(dat$start_position)
        mg=merge(dat,clinvar,by.x=xcols,by.y=ycols,all.x=T)
        
        matches=nrow(subset(mg,!is.na(variantid)))
        showNotification(sprintf("%s Clinvar matches",matches))
        
        return(mg)
      }
  }
}

# Function for getting website.
getWebsiteLink <- function(name){
  library(rvest)
  
  url = URLencode(paste0("https://www.google.com/search?q=",name))
  page <- read_html(url)
  results <- page %>% 
    html_nodes("cite") %>% # Get all notes of type cite. You can change this to grab other node types.
    html_text()
  result <- results[1]
  return(as.character(result)) # Return results if you want to see them all.
}


getCivicVariants<-function(){
  
  if (!file.exists("civic_vs.rds")){
  vs=fread("https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv")
  saveRDS(vs,"civic_vs.rds")
  }else{
    vs=readRDS("civic_vs.rds")
  }
  return(vs)
  
}


#Map data table to civic
mapToCivic<-function(dat="",vs=""){
  
  if (nrow(dat) >0 & nrow(vs)>0){
      vs=as.data.table(vs)
      vs$summary=NULL
    dat=as.data.table(dat)
    
    #print(names(dat))
    
    if ("chrom" %in% names(dat))
        dat$chrom=NULL
    if ("type" %in% names(dat))
       dat$type=NULL
  
    
    vs=merge(vs,dat,by.x=c("gene","variant"),by.y=c("hugo_gene_symbol","amino_acid_change_short"))
    vs=unique(vs,by=c("gene","variant","sample_name"))
    cols=names(vs)[1:13]
    cols=c("sample_name",cols)
    if ("coverage" %in% names(vs))
      cols=c(cols,"coverage")
    if ("allele_fraction" %in% names(vs))
      cols=c(cols,"allele_fraction")
 
    
  return(subset(vs,select=cols))
  }
}

variantTableUI<-function(id) {
  ns <- NS(id)
  uiOutput(ns("ui"))
}

#Test Application to view data in ES-------

variantTableTest<-function(dataobject=NULL){
  library(shinydashboard)
  smdf=as.data.table(getSampleIds())[order(-doc_count)]
  
  sampleids=smdf$sample_id
  names(sampleids)=paste0(smdf$sample_id,"(",smdf$doc_count,")")
  
  app <- shinyApp(
    ui=dashboardPage(skin = "black",title = "ES Variant Store Dashboard",
      dashboardHeader(title ="ES Variant Store Dashboard" ),
      dashboardSidebar(
           sidebarMenu(
             # Setting id makes input$tabs give the tabName of currently-selected tab
             id = "tabs",
             selectInput("qsample","Choose Sample",multiple = F,sampleids),
             menuItem("Variants", tabName = "dashboard", icon = icon("dashboard")),
             menuItem("Samples", tabName = "samples", icon = icon("dashboard")),
             menuItem("Deletions", tabName = "deletions", icon = icon("th")),
             menuItem("Data Status", icon = icon("th"), tabName = "datastat", badgeLabel = "new",badgeColor = "green")
           
           )
      ),
      dashboardBody(
        tabItems(
          tabItem("dashboard",
            variantTableUI("tb")
          )
          ,tabItem("samples",
                  DTOutput("smgenes")
          )
          ,tabItem("datastat",
                   DTOutput("filesloaded")
          )
          ,tabItem("deletions",
                   DTOutput("delstable")
          )
        )
      )
    ),
    server = function(input, output,session) {
      
      output$delstable<-renderDT({
      
       # if (!is.null(input$qsample)){
       #   dt=getDeletions(sid=input$qsample)
      #  }else{
        #  dt=getDeletions()
       # }
        dt=getDeletions()
        dt=as.data.table(dt)
        dt=dt[order(-reads)]
        datatable(dt,filter="top",rownames=F,options=list(scrollX=T,scrollY="600px"))
      })
      
      output$smgenes<-renderDT({
        req(input$qsample)
        dt=getSampleGeneCounts()
        dt=as.data.table(dt)
        datatable(dt,filter="top",options=list(scrollX=T,scrollY="600px"))
      })
      
      output$filesloaded<-renderDT({
        dt=getLoadedVariantFileStatus()
        dt=as.data.table(dt)
        dt$id=NULL
        datatable(dt,filter="top",options=list(scrollX=T,scrollY="600px"))
      })
      
      #Get data frame of variatns from elasticsearch
      dataobject<-reactive({
        req(input$qsample)
        showNotification("Getting data from ES...")
        t0=Sys.time()
        vdf=getVariants(qtype = "sample",qsample = input$qsample)
        t1=Sys.time()
        tres=difftime(t1,t0)
        showNotification(sprintf("%.2f sec Data ready,displaying",tres))
        printf("Query time for %s %.2f sec",input$qsample,tres)
        getSNVobject(vdf)
      })
      observe({
        input$qsample
        callModule(variantTable,"tb",dataobject(),genes=NULL,selcols=c("sample_name","amino_acid_change_short","location","ref","allele","hugo_gene_symbol","Variant","classification","coverage","allele_fraction","External Database Matches"))
    
      })  
        
      })
  runApp(app)
}


variantTable<-function(input,output,session,object,prefs=NULL,usePanelFilters=T,DomapToClinvar=FALSE,selcols=c("sample_name","hugo_gene_symbol","location","reference_allele","allele","amino_acid_change_short","Variant","sample_name","classification","External Database Matches"),genes=NULL,mycss=NULL,MAXLINES=500,...) {
  rv=reactiveValues(custompass=NULL,panelmatches=0)
  #Button listeneres
  
  #Append Custom Data Modal -----
  observeEvent(input$appendbtn,{
    ns=session$ns
    #dat =as.data.frame(alldata()) #get column names from reactive
    showModal(
      modalDialog(title = "Append Custom Database",fade = T,easyClose = T,size="m",  xbutton(),
                  shinyjs::useShinyjs(),
                  p("Upload data from a custom MS Excel sheet (.xls/.xlsx format only)."),
                  p("Variants will be joined based on matching data. Minimum information required is cNomen,Gene,Chrom,Transcript"),
                  fileInput(ns("customxl"),"Upload file",multiple=F, 
                            accept = c(".xlsx","application/vnd.ms-excel",
                                  "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                  ),
                  p("Please click \"Append\" to join the variant table data with the uploaded sheet."),
                  actionButton(ns("joindata"),"Append",icon=icon("play")),
                  actionButton(ns("resetjoin"),"Reset",icon=icon("reset")),
                  htmlOutput(ns("validatesheet")),
                  checkboxInput(ns("alljoin"),label="Show unmatched rows of data",value = T),
                  fluidRow(
                    column(9,
                      DTOutput(ns("appendsheet"))
                  )
                  ),
                  uiOutput(ns("appendedview"))
      )
    )
    shinyjs::disable("joindata")
    shinyjs::disable("resetjoin")
    jqui_resizable(selector = '.modal-content')
  })
  #Show error is sheet is invalid
  output$validatesheet<-renderText({
    invalidateLater(session = session,millis=3000)
    xldat=customdata()
  #  names(xldat)=gsub("\\s+","_",names(xldat))
    joinrequiredcols= c("Chrom","Gene","cNomen","Transcript")
    coltest=all(joinrequiredcols %in% names(xldat))
    ret=""
    if (coltest==TRUE){
      shinyjs::enable("joindata")
      shinyjs::enable("resetjoin")
      rv$custompass<-1
      ret="<i  style='color:blue' class=\"fa fa-3x fa-check-circle\"></i><h5 style='color:blue;font-weight: bold;'>Columns in input file are valid</h4>"
    }else{
      shinyjs::disable("joindata")
      shinyjs::disable("resetjoin")
      rv$custompass<-NULL
      ret="<i  style='color:red' class=\"fa fa-3x fa-times\"></i><h3 style='color:red;font-weight: bold;'>Columns not valid, need Chrom,Gene,cNomen and Transcript</h4>"
    }
    ret
      
  })
  observeEvent(input$resetjoin, {
    rv$custompass<-NULL
  })
  
  #get the gene panel data
  paneldata<-reactive({
    getPanel()
  })
  
  #turn the uploaded excel sheet into a data frame and modify some columns if they exist
  customdata<-reactive({
    req(input$customxl)
    if (grepl("xls$",input$customxl$name)){
      dat=readxl::read_xls(input$customxl$datapath)
      
    }else if (grepl("xlsx$",input$customxl$name)){
      dat=readxl::read_xlsx(input$customxl$datapath)
    }
    #dat=as.data.table(dat)
    if ("Chrom" %in% names(dat)){
     dat$Chrom=gsub("chr","",dat$Chrom)
      dat$Chrom=paste0("chr",dat$Chrom)
    }
    if ("SEQ" %in% names(dat)){
    dat$REF=gsub("REF=","",str_extract(dat$SEQ,"REF=[^;]+"))
    dat$ALT=gsub("OBS=","",str_extract(dat$SEQ,"OBS=[^;]+"))
    }
    return(dat)
  })
  
  #preview of the table
  output$appendsheet<-renderDT({
     req(input$customxl)
    
     dat=customdata()
     if (nrow(dat)>100){
       dat=head(dat,100)
     }
     if (!is.null(rv$custompass))
       datatable(dat,caption="Uploaded file preview",filter="top",options=list(scrollX=T,autowidth=T,pageLength=5))
  })
  
  
  
  #Custom Column selector modal ------
  observeEvent(input$colbtn,{
    ns=session$ns
    dat =as.data.frame(alldata()) #get column names from reactive
    showModal(
      modalDialog(title = "Choose Columns to Display",fade = T,easyClose = T,size="m",  xbutton(),
                  textInput(ns("colsearch"),"Search",value=NULL),
                uiOutput(ns("dispcols"))
      )
    )
    jqui_resizable(selector = '.modal-content')
  })
  
  output$dispcols<-renderUI({
    ns=session$ns
   # req(input$colsearch)
    dat =as.data.frame(alldata()) #get column names from reactive
    nms=names(dat)
    if (is.null(input$colsearch)){
      nms=names(dat)
    }else{
    nms=grep(input$colsearch,nms,value=T)
    }
    checkboxGroupInput(ns("displaycols"),
                       h5("Show Table Columns:"),
                       selected=selcols,
                       sort(as.character(nms)),
                       width="70%",inline=T
    )      
  })
  
  #Variant exclusion modal -------
  observeEvent(input$exlbtn,{
    ns=session$ns
    dat =as.data.frame(alldata()) #get column names from reactive
    varclasses=unique(dat$classification)
    fvarclasses=grep("^synonymous|intron|UTR|upstream|downstream",varclasses,value=T)
    showModal(
      modalDialog(title = "Variant Class Exclusion Filters",fade = FALSE,easyClose = T,size="m",  xbutton(),
                  checkboxGroupInput(ns("varclass"),
                                     h5("Exclude by Variant Class:"),
                                     sort(as.character(varclasses)),
                                     selected=fvarclasses,
                                     width="70%",inline=T
                                     
                  )
                  )       
      )
    jqui_resizable(selector = '.modal-content')
  })
  
  #Modal for Preset filters --------
  observeEvent(input$presetbtn,{
    ns=session$ns
    pn=paneldata()
    panels=panelNames(paneldata())
    dat =as.data.frame(alldata()) #get column names from reactive
    showModal(
      modalDialog(title = "Preset Filters",fade = T,easyClose = T,size="m",  xbutton(),
                  fluidPage(
                    selectInput(ns("panelid"),"Choose Panel",multiple=F,c("None",panels)),
                    actionButton(ns("usepanel"),"Go",icon=icon("play"))                  
                    )
      )
    )
    jqui_resizable(selector = '.modal-content')
  })
  
  #Modal for something else not used -------
  observeEvent(input$advbtn,{
    ns=session$ns

    dat =as.data.frame(alldata()) #get column names from reactive
    showModal(
      modalDialog(title = "Advanced Panel Filters",fade = T,easyClose = T,size="m",  xbutton(),
            "Not enabled yet"
            
      )
    )
    jqui_resizable(selector = '.modal-content')
  })
  
if (!is.null(mycss))
  mycss="border: 1px solid silver;overflow-y:scroll; overflow-x:scroll ; tbody {  padding: 5px 2px 2px 1px; font-family: tahoma; font-weight:normal; font-size:x-small}"
  #cstyle='overflow-x:scroll;height:800px;max-height:"1200px";font-size:inherit;border: 1px solid silver;' 
###Render the UI
  
output$activebutton<-renderUI({
  ns=session$ns
  pvar=NULL
  if (!is.null(input$panelid)){
    pvar=input$panelid
    col="blue"
    if (input$panelid=="None"){
      pvar=NULL
    }
  }
  if (!is.null(pvar) && input$usepanel>0){
    div(title="Active Gene Panel",actionButton(ns("unusdbtsrsdsds"),pvar,class="btn btn-lg",style="color:white;background-color:orange;font-weight:bold;"))
  }

})
  
 #Main Data table for our data that's dynamically filtered -----
output$ui <- renderUI({
  ns=session$ns
  dat =as.data.frame(alldata()) #get column names from reactive
  varclasses=unique(dat$classification)
  fvarclasses=grep("^synonymous|intron|UTR|upstream|downstream|non_coding",varclasses,value=T)
  
 #Start with some buttons
  menubuttons=div(id=ns("fileactionbuttons"),class="btn-group float-right",role="group",
   actionBttn(ns('colbtn'),size="xs",style="pill",label=div(title="Choose Columns","Choose Columns"),icon=icon("plus")),
   actionBttn(ns('appendbtn'),size="xs",style="pill",label=div(title="Append Custom Database","Custom Database"),icon=icon("plus")),
   
   actionBttn(ns('exlbtn'),size="xs",style="pill",label=div(title="Exclusion Filters","Exclusion Filters"),icon=icon("filter")),
   actionBttn(ns('presetbtn'),size="xs",style="pill",label=div(title="Preset Rules","Preset Rules"),icon =icon("list")),
   actionBttn(ns('advbtn'),size="xs",style="pill",label=div(title="Advanced Filters","Advanced Filters"),icon=icon("cog"))
  )
  
  fluidPage(
    fluidRow(
      menubuttons
    ),
  ############# Filters #####
    fluidRow(
     column(1,selectInput(ns("genefilter"),"Gene Filter",unique(dat$hugo_gene_symbol),multiple=T,selected=NULL )),
     column(1,selectInput(ns("samplefilter"),"Sample Filter",gsub(".vcf|.tsv|.info|.gz","",object$getSampleIds()),selected=NULL,multiple=T )),
     column(1,textInput(ns("hgvspattern"),label="HGVSc/HGVSP/Loc.",value=NULL)),
     column(1,selectInput(ns("snvtype"),label="Variant Type",multiple=T,unique(dat$type),selected=NULL )),
     column(1, uiOutput(ns("activebutton"))),
     column(1,checkboxInput(ns("showonc"),label="OncoKB Matches only",value = FALSE)),
     column(1,checkboxInput(ns("uniqvars"),label="Unique Variants only",value = FALSE)),
     column(1,materialSwitch(inputId = ns("filtersoff"), label = "Deactivate class filters",value=F, status = "danger")),
     #Stats
     column(3, 
            div(class="span  label-warning",h4(countupOutput(ns("varcount")))),
            div(class="span  label-info",h4(countupOutput(ns("recordcount")))),
            countupOutput(ns("samplecount")),"/",countupOutput(ns("totalsamplecount"))
            )
     ),
  fluidRow(
    tabsetPanel(
    tabPanel(h5("Table"),  
    downloadLink(ns("downloadtable"),'Download CSV'),
    column(12,DT::dataTableOutput(ns("table")),style=mycss)
    )
    ,tabPanel(h6("OncoKB"),
              div(DT::dataTableOutput(ns("onctable")),style=mycss)
    ),
    tabPanel(h6("CIViC"),
             div(DT::dataTableOutput(ns("civictable")),style=mycss)
    )
    )
  )
  )
})

output$civictable<-renderDT({
  tbl<-civicdat()
  tbl=as.data.table(tbl)
  #tbl$l
  tbl[,variant:=sprintf('<a target="_blank" href="%s">%s</a>',variant_civic_url,variant)]
  #setnames(tbl,"l","variant")
  tbl$variant_civic_url=NULL
  
  datatable(tbl,caption="matches to Civic",filter="top",escape = F,options=list(
    scrollX=T,
    dom="pBtrilf",
    pageLength=nrow(tbl)
  )
  )
})

output$onctable<-renderDT({
  dat=object$getTable(doFilter=FALSE)
  dat=as.data.table(dat)
  if ("oncokb" %in% names(dat)){
    dat<-subset(dat,!is.na(oncogenicity))
    print(nrow(dat))
    dat$sample_id=gsub(".fastq|.vardict|.qc|.gz|.bam","",dat$sample_id)
    cols=c("hugo_gene_symbol","amino_acid_change_short","isoform","oncogenicity","mutation_effect")
    summ<-dat[hugo_gene_symbol!="-",by=cols,list(
      count=length(unique(sample_id)),
      Samples=paste0(unique(sample_id),collapse=",")
    )][order(-count)]
    setnames(summ,"count","No Of Samples")
    datatable(summ,caption="matches to OncoKB",filter="top",options=list(
      scrollX=T,
      dom="pBtrilf",
      pageLength=nrow(summ)
    )
    )
  }else {
    summ=data.frame(Result="No Matching Records")
    datatable(summ)
  }

})


output$totalsamplecount<-renderCountup({
  opts = list(
    useGrouping = TRUE, 
    separator = ',',
    suffix=" Sample(s)"
  )
  
  dat=alldata()
  cnt=0
  if ("sample_name" %in% names(dat)){
    cnt=length(unique(dat$sample_name))
  }

  countup(count = cnt, start = 0,    
          height='200px',
          width='200px',
          duration=0.8, 
          options = opts)
  
})

output$samplecount<-renderCountup({
  dat=allVarTable()
  dat=subset(dat,!is.na(sample_name))
  val=length(unique(dat$sample_name))
  
  opts = list(
    useEasing=T,
    useGrouping = TRUE, 
    separator = ',',
    suffix=""
  )
  countup(count = val, start = 0,    
          height='200px',
          width='200px',
          duration=0.8, 
          options = opts)
  
})

output$recordcount<-renderCountup({
  opts = list(
    useEasing=TRUE,
    useGrouping = TRUE, 
    separator = ',',
    suffix=" Records"
  )
  countup(count =nrow(allVarTable()) , start = 0,    
          height='200px',
          width='200px',
          duration=0.8, 
          options = opts)
  
})

output$varcount<-renderCountup({
  opts = list(
    useEasing=TRUE,
    useGrouping = TRUE, 
    separator = ',',
    suffix=" Unique Variants"
  )
  countup(count =length(unique(allVarTable()$location)) , start = 0,    
          height='200px',
          width='200px',
          duration=0.8, 
          options = opts)
  
})

output$selrows<-renderPrint({
    tb=allVarTable()
    if (!is.null(input$table_rows_selected))
    subset(tb[input$table_rows_selected,],select=c("hugo_gene_symbol","sample_name","variant_aa"))
})

#Add colored annotatons  to main variant Table
addAnnotations<-function(dat) {
  #Add these new columns
  dat$Annotation=""
  dat$adbsnp=""
  dat$acos=""
  dat$ads=""
  dat$akg=""
  dat$amills=""
  dat$aexac=""
  #Markup dbsnp annoation
  if ("is_dbsnp" %in% names(dat) ) {
    dbsnpurl="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs"
    if ("external_annotation_dbsnp" %in% names(dat) ) {
      dat[grepl("^1$",is_dbsnp,perl=T,ignore.case=T),
          adbsnp:=paste(sprintf('<a target="_blank" type="submit" class="btn btn-primary btn-xs" href="%s=%s">dbSNP</a>',
                                dbsnpurl,external_annotation_dbsnp),sep="")]
    }else {
      if ("VariationID" %in% names(dat)){
       dat[grepl("^1$",is_dbsnp,perl=T,ignore.case=T),
          adbsnp:=paste(adbsnp,
                        sprintf('<a target="_blank" href="https://www.ncbi.nlm.nih.gov/clinvar/variation/%s" type="submit" type="button" class="btn btn-primary btn-xs" >dbSNP</a>',VariationID),
                                sep=" ")]
      }else{
        dat[grepl("^1$",is_dbsnp,perl=T,ignore.case=T),
            adbsnp:=paste(adbsnp,
                          sprintf('<a target="_blank" href="https://www.ncbi.nlm.nih.gov/clinvar/variation/%s" type="submit" type="button" class="btn btn-primary btn-xs" >dbSNP</a>',"NA"),
                          sep=" ")]
        
      }
      
    }
    dat$adbsnp=gsub("^rs","",dat$adbsnp)	
  }
  
  if ("is_mills" %in% names(dat)) {
    dat[grepl("1",is_mills),amills:=paste(amills,' <button type="button" class="btn btn-warning btn-xs">Mills</button>',sep=" ")]
    
  }
  
  if ("is_1kgsnp" %in% names(dat)) {
    dat[grepl("1",dat$is_1kgsnp,perl=T,ignore.case=T),akg:=paste(akg,' <button type="button" class="btn btn-warning btn-xs">1KG</button>',sep=" ")]
  }
  
  if ("is_exac" %in% names(dat)) {
    dat$exacdat=paste(dat$chrom,dat$pos,dat$reference_allele,dat$variant_allele,sep="-") #Make new coluhmn  on the fly
    dat[grepl("1",is_exac,perl=T,ignore.case=T),
        aexac:=paste(sprintf('<a target="_blank" type="submit" class="btn btn-warning btn-xs" href="http://exac.broadinstitute.org/variant/%s">ExAc</a>',exacdat),sep="")]
    
  }
  
  if ("dbsnp_somatic" %in% names(dat)) {
    dat[grepl("1",dbsnp_somatic,perl=T,ignore.case=T),ads:=paste(ads,' <button type="button" class="btn btn-warning btn-xs">dbSNP-somatic</button>',sep=" ")]
    
  }
  
  #TODO fix this markup problem
  if ("is_cosmic" %in% names(dat) & "cosmicidstr" %in% names(dat)) {
    dat[grepl("1",is_cosmic,perl=T,ignore.case=T),acos:=paste(
      sprintf('<a target="_blank" type="submit" class="btn btn-danger btn-xs" href="http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=%s">COSMIC</a>',
              sub("COSM","",cosmicidstr)
      ),sep="")]
  }
  
  if ("oncokb" %in% names(dat)) {
    dat$adonc=""
    dat[oncokb==TRUE,adonc:=paste(adonc,sprintf('<a type="button" target="_blank" type="submit" class="btn btn-warning btn-xs" href="http://oncokb.org/#/gene/%s/alteration/%s">OncoKB</a>',symbol,amino_acid_change_short),sep=" ")]
  }
  
  #onclick = 'Shiny.onInputChange(\"killjob_btn\",  this.id); Shiny.onInputChange("killjob_click", Math.random())'
  
  #dat$Annotation=paste(dat$adbsnp,dat$acos,dat$aexac,sep=" ")
  dat$Annotation=paste(dat$adbsnp,dat$acos,dat$ads,dat$akg,dat$amills,dat$adonc,dat$aexac,sep=" ")
  return(dat)
}


addURLs <- function(dat) {
  dat$amino_acid_change=dat$variant_aa
  exacu="http://exac.broadinstitute.org/variant/"  
  dat$exacdat=paste(dat$chrom,dat$pos,dat$reference_allele,dat$variant_allele,sep="-")
  dat$exacurl=sprintf("<a href=\"%s\" target=\"_blank\" >%s</a>",paste(exacu,dat$exacdat,sep=""),dat$variant_aa)
  dat$variant_aa=dat$exacurl
  dbu="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs="
  
  
  if ("external_annotation_dbsnp" %in% names(dat) ) {
    dat$rsid=sub("rs","",dat$external_annotation_dbsnp)
    dat$rsurl=sprintf("<a href=\"%s\" target=\"_blank\" >%s</a>",paste(dbu,dat$rsid,sep=""),dat$external_annotation_dbsnp)
    dat$variant_aa= paste(dat$variant_aa,dat$rsurl,sep="<br>")
  }
  
  if ( "cosmicidstr" %in% names(dat))  {
    cosurl="http://cancer.sanger.ac.uk/cosmic/mutation/overview?id="
    dat$cid=dat$cosmicidstr
    dat$cid=gsub(",COSM\\S+","",dat$cid,perl=T)
    dat$cid=sub("COSM","",dat$cid,perl=T)
    dat$cz=dat$cid
    dat$cid=sprintf("<a href=\"%s\" target=\"_blank\" >COSM%s</a> ",paste(cosurl,dat$cid,sep=""),dat$cid)
    #dat[dat$cz!="0",][,"variant_aa"] = paste(dat$variant_aa,dat$cid,sep="<br>")
    dat$variant_aa= paste(dat$variant_aa,dat$cid,sep="<br>")
  }
  
  return(dat)
  
}


variantclickback=function(){JS('function(table) {
                    table.on("click.dt","tr", function() {
                    var data=table.row(this).data();    
                    console.log(jQuery(data[1]).text() + " Gene Clicked on gene , Variant is:" + data[2] + " Sample " + data[0]);
                    Shiny.onInputChange("variantclickgene",jQuery(data[1]).text());
                    Shiny.onInputChange("variantclickaa",data[2]);	
                    Shiny.onInputChange("variantclicksm",data[0]);
                    });}')
}

#All data (no filtering) -----
alldata<-reactive({
  dat=object$getTable(doFilter=FALSE)
  dat$variantid=dat$variant_aa
  dat=as.data.table(dat)
  setnames(dat,tolower(names(dat)))
  dat
})

# Civic data mapping --------
civicdat<-reactive({
  mapToCivic(dat=alldata(),vs=getCivicVariants())
})

#Reactive data with filters -------
allVarTable<-reactive({
  dat=object$getTable(doFilter=FALSE)
  dat$variantid=dat$variant_aa
  dat=as.data.table(dat)
  #Add some columns
  if ("CHROM" %in% names(dat)) {
    dat$CHROM=NULL
  }
  if ("TYPE" %in% names(dat)){
    dat$TYPE=NULL
  }
  #print(names(dat))  
  setnames(dat,tolower(names(dat)))
  
  # saveRDS(dat,"vtabledata.rds")
  
   #Uniquify variants
  showunique=TRUE
  
  #Show Unique variants only
  if (!is.null(input$uniqvars)) {
      if (input$uniqvars ==TRUE){
         ucols=c("sample_name","hugo_gene_symbol","amino_acid_change_short","location")
        dat=unique(dat,by=ucols)   
    }
  }


  if (!is.null(input$filtersoff)){
    if (isTRUE(input$filtersoff)){
      showNotification("Not doing any filtering...",type="message")
    }else{
      if (!is.null(input$varclass)) {
        dat=subset(dat,!(classification %in% input$varclass) )
      }else{
        dat=subset(dat,!(grepl("^synonymous|intron|UTR|upstream|downstream",classification)))
        showNotification("Applying default upstream|intron|downstream|synonymous|UTR filters",type="warning")
      }
    }
  }
  


  if (!is.null(input$genefilter)) {
   # if (!is.na(input$genefilter))
     dat=subset(dat,hugo_gene_symbol %in% input$genefilter )
     printf("%s rows returned from gene filter search",nrow(dat))
  }
  
  if (!is.null(input$samplefilter)) {
    if (!is.na(input$samplefilter))
      dat=subset(dat,sample_name %in% input$samplefilter )
  }
  
  if (!is.null(input$snvtype)) {
    if (!is.na(input$snvtype))
      dat=subset(dat,type %in% input$snvtype )
  }
  
  #General purpose filter
  if (!is.null(input$hgvspattern)){
    if (input$hgvspattern !="")
      dat=subset(dat,grepl(input$hgvspattern,hgvsc) | grepl(input$hgvspattern,hgvsp) |  grepl(input$hgvspattern,location) | grepl(input$hgvspattern,amino_acid_change_short) )
  }
  
  #Show ONcoKB matches
  if (!is.null(input$showonc)) {
    if ("oncokb" %in% names(dat)) {
      if (input$showonc==TRUE)
         dat=subset(dat,oncokb==input$showonc)
    }
  }
  
  
  
  requiredcols=c("hugo_gene_symbol","variant_aa","sample_name","label")
  
  if (! unique((requiredcols %in% names(dat)))   ) {
       showNotification(sprintf("Input data failed on column check, need %s",paste0(requiredcols,collapse="\n")),duration=6,type="error")
       return(NULL)
  }
  
  if (!is.null(genes)) { 
      dat=subset(dat,hugo_gene_symbol %in% genes)
  }
  
  dat$url=sprintf("<a href=\"../samples/s=%s\">%s</a>",dat$sm,dat$sm)
  dat$geneurl=sprintf("<a href=\"http://www.ncbi.nlm.nih.gov/gene/?term=%s[sym]\" target=\"_blank\" >%s</a>",dat$gene,dat$gene)
  dat$gene=dat$geneurl
  dat$Annotation=""
  
  
    if (DomapToClinvar==TRUE){
      showNotification("Mapping to clinvar")
      dat=mapToClinvar(dat=dat)
    }
  

  tblcopy=copy(dat)
  
  #Preset gene panel filtering
  if (usePanelFilters==TRUE){
    if (!is.null(input$usepanel)){
      if (input$usepanel >0){
        pnid=isolate(input$panelid)
        if (pnid !="None"){
          print("Applying gene panel filters")
          printf("Panel=%s",pnid)
          gdat=filterByPanel(indata=dat,paneldata=paneldata(),panelname = pnid)
          showNotification(sprintf("%s result rows returned from Panel %s",nrow(gdat),pnid),type="warning")
          #replace the data
          rv$panelmatches=nrow(gdat)
          rv$currpanel=pnid
          dat=gdat
        }
      }
    }
  }
  
  
  
  
  #Join Custom Database actions
  if (!is.null(input$joindata)){
    showall=F
    if (!is.null(input$alljoin)){
      if (isTRUE(input$alljoin)){
        showall=T
      }else{
        showall=F
      }
    }
    
    if (input$joindata >0){
   # req(input$customxl)
    showNotification("Custom Data Found, appending",type="warning",duration=4)
    xldat=customdata()
    names(xldat)=gsub("\\s+","_",names(xldat))
    joinrequiredcols= c("Chrom","Gene","cNomen","Transcript")
    coltest=all(joinrequiredcols %in% names(xldat))
    if (coltest==FALSE){
      showNotification("Insufficient columns to join data",type="error")
    }
    if (coltest==TRUE){
         dat= joinCustomDb(dat=dat,xldat=xldat,showall=showall)
    }
    }
  }
  #end of custom data join
  if (!is.null(input$resetjoin)){
    if (input$resetjoin>0)
        dat=tblcopy
  }
  
  ## Add annotations
  if (1==1) {
    dat=addAnnotations(dat)  #add button annotations
    dat=addURLs(dat)
  }
  
  
  dat$variant_aa=gsub("COSM0","",dat$variant_aa)
  dat$variant_aa=gsub("\\>\\.","",dat$variant_aa,perl=T)
  
  if ("gene" %in% names(dat) && !("Gene" %in% names(dat))){
      setnames(dat,"gene","Gene")
    dat<-dat[order(Gene)]
  }
  
  dat$Variant=dat$variant_aa
  setnames(dat,"Annotation","External Database Matches")
  
  dat[order(-hugo_gene_symbol)]
  
})


#Download the table --------
output$downloadtable<- downloadHandler(
  filename= function(){
       paste('variantdata-', Sys.Date(), '.csv', sep='')
  },
  content=function(file){
      mytable= allVarTable()
    
      staticols=c('variantid',selcols)
      if (!is.null(input$displaycols)) {
        mytable=subset(mytable,select=unique(c(staticols,input$displaycols)))
      }
      mytable[["External Database Matches"]]=NULL
      mytable$Variant=NULL
      write.csv(mytable,file,row.names=F,quote=T)
    }
)
  
#The actual main table ID --------  
output$table <- DT::renderDataTable({
  ns=session$ns
  mytable= allVarTable()

  staticols=selcols
  
  #Add some custom columns
  if ("mutation effect" %in% names(mytable)){
    staticols=c(staticols,"mutation effect") 
  }
  
  #Add custom data columns
  if (!is.null(input$customxl)){
    xldat=customdata()
    matchingnames=intersect(names(xldat),names(mytable))
    print(matchingnames)
    staticols=c(staticols,matchingnames)
    MAXLINES=1000
  }
  
  #print(names(mytable))
  #Set the display columns to show in main view
  if (!is.null(input$displaycols)) {
    mytable=subset(mytable,select=unique(c(staticols,input$displaycols)))
  }else{
    mytable=subset(mytable,select=staticols)
  }
    #Remov the underscroes from table names
  setnames(mytable,gsub("_"," ",names(mytable)))

  mytable=head(mytable,MAXLINES)
  
  mycaption="Variants"
  pnid=isolate(input$panelid)
  if (!is.null(pnid)) {
    if (pnid!="None"){
      mycaption=paste0("Variants matching ",pnid)
    }
  }
  
  
  #Add View button
  Actions = paste(
    shinyInput(actionLink, nrow(mytable), 'button_', label = icon("info-circle"),
               onclick = sprintf('Shiny.onInputChange(\"%s\",  this.id);
                                 console.log("Var click" + String(this.id));
                                 Shiny.onInputChange(\"%s\",  Math.random());',
                                 ns("vdetail_btn"),
                                 ns("vdetail_click")),
               style='padding:0px;' )
  )
  #suppressWarnings({ Actions=gsub("submenow",icons,Actions) })
  mytable=as.data.table(mytable)
  mytable=cbind(Actions,mytable)
  setnames(mytable,"Actions",".")
  
  #Some formatting
  if ("classification" %in% names(mytable)){
    mytable$classification=gsub(",","\n",mytable$classification)
    mytable$classification=gsub("_"," ",mytable$classification)
  }
  if ("Variant" %in% names(mytable)){
    mytable$Variant=gsub(":","\n",mytable$Variant)
    mytable$Variant=gsub("COSMNA|\\s\\-\\s","",mytable$Variant)
    
  }
  
  cn = stringr::str_wrap(colnames(mytable),width = 10)
  
  #Add ellipsis truncation for ref/alt if they're too long
  truncCols = list(list(
    targets  = c(5,6),
    render   = JS(
      "function(data, type, row, meta) {",
      "return type === 'display' && data.length > 6 ?","'<span title=\"' + data + '\">' +
      data.substr(0, 6) + '...</span>' : data;", "}")))
  
  setnames(mytable,cn)
  
  datatable(mytable,escape=F,filter="bottom",rownames=T,caption=mycaption,
            extensions =  c('ColReorder','Buttons'),
    options = list(
      buttons = c('copy','colvis'),
    colReorder = TRUE,
    scrollX=TRUE,
    scrollY="450px",
    autoWidth=T,
    columnDefs =truncCols,
    pageLength = 100, 
    lengthMenu = c(15,20,50,100,200,nrow(mytable)),
    dom='pBTlifrt',
    initComplete = JS(
               "function(settings, json) {",
               "$(this.api().table().header()).css({'font-size': '79%','font-weight':'normal','background-color':'#3c8dbc','color':'white'});",
               "$(this.api().table().body()).css({'font-size': '85%','font-weight':'normal','background-color':'white'});",
               "}")
    )
  )  %>%  formatStyle(columns = names(mytable)[1] , width='1%')
  #%>% formatStyle(columns = names(mytable), fontSize = '83%')
  
})


#Modal for the variant detailed view
variantModal<-modalDialog(title="Variant Details",
                            xbutton(),
                            hr(),
                            fluidPage(
                            uiOutput(session$ns("vdetailpage")),
                            DTOutput(session$ns("vdetailtbl"))
                            ),
                            easyClose = T,fade=F,size="l"
)

#Variant details modal popup
output$vdetailpage<-renderUI({
  dat=rv$vardatarow
  dat=head(dat,1)
  readcount=NA
  af=NA
  #vdetail=paste(dat$location)
  if (!is.null(dat$allele_fraction)){
    af=dat$allele_fraction
  }
  if (!is.null(dat$coverage)){
    readcount=dat$coverage
  }
  tagList(
   fluidRow(
     shinydashboard::infoBox("Gene",value=ifelse(!is.null(dat$hugo_gene_symbol),dat$hugo_gene_symbol,""),icon=icon("cubes"),width=3),
     shinydashboard::infoBox("Sample",value=dat$sample_name,icon=icon("vial"),width=3),
    fluidRow(
    shinydashboard::valueBox(prettyNum(readcount,big.mark = ","),color="orange", "Reads", icon = icon("dna"),width=3)
    ,shinydashboard::valueBox(sprintf("%.3f",af),color = "purple", "Allele Freq", icon = icon("check"),width=3)
    ),
    fluidRow(
    shinydashboard::infoBox("Variant",value=paste(dat$ref,">",dat$allele,"\n",dat$amino_acid_change_short),icon=icon("info"),width=6),
    shinydashboard::infoBox("Position",value=paste(dat$location),icon=icon("info"),width=4,color="green")
    )
    )
   )
    
    
})

output$vdetailtbl<-renderDT({
  if (!is.null(rv$vardatarow)) {
    dat<-rv$vardatarow
    dat<-t(dat)
    dat=as.data.frame(dat)
    dat$Field=rownames(dat)
    names(dat)<-c("Value","Field")
    dat<-subset(dat,Value !="-" ,select=c("Field","Value"))
    dat=subset(dat,!is.na(Value) | Value !="")
    datatable(dat,filter="top",selection="none",escape=F,rownames=F,caption="Variant Details",options=list(scrollX=T,scrollY="400px",pageLength=nrow(dat)))
  }
  
})

#Click event on table View button
observeEvent(input$vdetail_click, {
  selectedRow <- as.numeric(strsplit(input$vdetail_btn, "_")[[1]][2])
  tbl=as.data.frame(allVarTable())
  datarow=tbl[selectedRow,]
  rv$vardatarow=datarow
 showModal(variantModal)
  jqui_resizable(selector = '.modal-content')
})


## END OF APP MODULE -----

}

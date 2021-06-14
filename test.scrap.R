
sampCounts <- list()
for( ss in samplist ){
  message("-",ss);
  sampCounts[[ss]] <- list();
  kk <- dda$sample == ss & dda$IS_SCNA == 1;
  ddax <- dda[kk,]
        kkn <-  dda$sample == ss & (dda$IS_SCNA == 1 | dda$event.type == "NA" | dda$event.type == "???");
        ddaxn <- dda[kkn,]
  multidata.info[[ p(ss) ]] <- list()
  multidata.info[[ p(ss) ]][["data"]] <- ddaxn;
  multidata.info[[ p(ss) ]][["events"]] <- character();
  multidata.info[[ p(ss) ]][["events.clonal"]] <- character();
  multidata.info[[ p(ss) ]][["events.info"]] <- list();
  multidata.info[[ p(ss) ]][["events.idx"]]  <- list();
  multidata.info[[ p(ss) ]][["events.data"]] <- list();
  multidata.info[[ p(ss) ]][["fullchrom"]] <- list();
  multidata.info[[ p(ss) ]][["flags"]] <- c();

    sampCounts[[ss]][[p("rawCT.ALL")]] <- sum(kk);
    sampCounts[[ss]][[p("chrCT.ALL")]] <- 0;
    sampCounts[[ss]][[p("mrgCT.ALL")]] <- 0;
    for(cc in chromlist){
      kkx  <- kk  & dda$Chromosome == cc
      kkxc <- kkx & dda$event.clonality == "clonal"
      if(sum(kkx) > 0){
        sampCounts[[ss]][[p("chrCT.ALL")]] <- sampCounts[[ss]][[p("chrCT.ALL")]]+1;
      }
    }
  #for(ee in eventlist ){
  #  kk <- dda$sample == ss & dda$event.type == ee & dda$IS_SCNA == 1;
  #  print(dda[kk,])
  #}

  for(ee in eventlist ){
    kk <- dda$sample == ss & dda$event.type == ee & dda$IS_SCNA == 1;

    sampCounts[[ss]][[p("rawCT.",ee)]] <- sum(kk);
    sampCounts[[ss]][[p("chrCT.",ee)]] <- 0;
    sampCounts[[ss]][[p("mrgCT.",ee)]] <- 0;
    #sampCounts[[ss]][[p("rawCT.clonal.",ee)]] <- sum(kk);
    #sampCounts[[ss]][[p("chrCT.clonal.",ee)]] <- 0;
    sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- 0;
    for(cc in chromlist){
      kkx  <- kk  & dda$Chromosome == cc
      kkxc <- kkx & dda$event.clonality == "clonal"
      if(sum(kkx) > 0){
        sampCounts[[ss]][[p("chrCT.",ee)]] <- sampCounts[[ss]][[p("chrCT.",ee)]]+1;
        idx <- which(kkx | (dda$sample == ss & (dda$event.info == "NA" | dda$event.type == "???") & dda$Chromosome == cc))
        while( length(idx) > 0 && ("NA" == ( dda$event.type[idx[1]] ) | "???" == ( dda$event.type[idx[1]] ) )     ){
          idx <- idx[-1];
        }
        while( length(idx) > 0 && ("NA" == ( dda$event.type[[ idx[[length(idx)]] ]] ) | "???" == ( dda$event.type[[ idx[[length(idx)]] ]] )  )){
          idx <- idx[-length(idx)];
        }


        if(length(idx) == 1){
          sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]+1;
          if( sum(kkxc) > 0){
            sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- sampCounts[[ss]][[p("mrgCT.clonal.",ee)]]+1;
            multidata.info[[ p(ss) ]][["events.clonal"]] <- c(multidata.info[[ p(ss) ]][["events.clonal"]],p(ee,":chr",cc))

          }
          is.clonal <- sum(kkxc) > 0
          clonality.string <- ifelse(sum(kkxc) > 0,"clonal","nonclonal")

            ddax <- multidata.info[[ p(ss) ]][["data"]]
            zzix <- which( ddax$rownum == dda$rownum[ idx[1] ])
            ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("SINGLETON[",ee,"]"))
            multidata.info[[ p(ss) ]][["data"]] <- ddax;
            multidata.info[[ p(ss) ]][["events"]] <- c(multidata.info[[ p(ss) ]][["events"]],p(ee,":chr",cc))
                  multidata.info[[ p(ss) ]][["events.info"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- list(
                      event = ee,chrom=cc,
                      start = ddax$Start_bp[[ zzix ]],
                      end   = ddax$End_bp[[zzix]],
                      clonal = is.clonal,
                      idx    = p(zzix,collapse=","),
                      N = 1
                  )
                  #message("class: ",class( multidata.info[[ p(ss) ]][["events.idx"]] ))
                  multidata.info[[ p(ss) ]][["events.idx"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- c(zzix);
                  multidata.info[[ p(ss) ]][["events.data"]][[ length(multidata.info[[ p(ss) ]][["events.data"]])+1  ]] <- ddax[zzix,,drop=F]
        } else if(length(idx) > 1){
          ddax <- multidata.info[[ p(ss) ]][["data"]]
          message("Multi On Chrom: ",ss,"/",cc);
          multi.on.chrom.ct = multi.on.chrom.ct  + 1;
          sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]+1;
          zzix <- which( ddax$rownum == dda$rownum[ idx[1] ])
          ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("START[",ee,"]"))
          if(DEBUG.MODE){ message("notes[i=",1,"]/[idx=",idx[1],"]/[ddax=",(zzix),"]=END") }
          multidata.info[[ p(ss) ]][["events"]] <- c(multidata.info[[ p(ss) ]][["events"]],p(ee,":chr",cc))
          #multidata.info[[ p(ss) ]][["events.info"]] <- c(multidata.info[[ p(ss) ]][["events.info"]],p(ee,":chr",cc,":", ddax$Start_bp[[zzix]],":"  ))

          for(i in seq_along(idx)[-1]){
            if(DEBUG.MODE){ message("i=",i) }
            if(DEBUG.MODE){ print(dda[idx[(i-1):(i)],c("sample","Chromosome","Start_bp","End_bp","event.type","rownum","notes")] )}
            if( dda$End_bp[idx[i-1]] + BUFFER.WINDOW < dda$Start_bp[idx[i]] || dda$rownum[idx[i-1]] + 1 < dda$rownum[idx[i]] ){
              if(DEBUG.MODE){ message("BREAK-step 1") }
              zzix <- which( ddax$rownum == dda$rownum[ idx[i-1] ])
              ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("END[",ee,"]"))
              if(DEBUG.MODE){ message("notes[i=",i-1,"]/[idx=",idx[i-1],"]/[ddax=",(zzix),"]=END") }

              idx.prev   <- 1:(i-1);
              zzix.prev  <- which( ddax$rownum %in% dda$rownum[ idx[idx.prev] ] )
              notes.prev <- ddax$notes[ zzix.prev ];
              ixx.start  <- max( which(grepl(p("START[",ee,"]"),notes.prev,fixed=T) ))
              idx.curr   <- ixx.start:(i-1);
              zzix.curr  <- which( ddax$rownum %in% dda$rownum[ idx[idx.curr] ])
              zzix.last  <- zzix.curr[[length(zzix.curr)]]
              is.clonal  <- any( f.na( ddax$event.clonality[zzix.curr] == "clonal" ))
              is.event   <- any( f.na( ddax$event.type[zzix.curr] == ee ))
              clonality.string <- ifelse(is.clonal,"clonal","nonclonal")
              if(is.clonal & is.event){
                sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- sampCounts[[ss]][[p("mrgCT.clonal.",ee)]]+1;
                multidata.info[[ p(ss) ]][["events.clonal"]] <- c(multidata.info[[ p(ss) ]][["events.clonal"]],p(ee,":chr",cc))
              }
              if(! is.event ){
                message("DELETING EVENT!")
                sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]-1;
                multidata.info[[ p(ss) ]][["events"]] <- multidata.info[[ p(ss) ]][["events"]][- length(multidata.info[[ p(ss) ]][["events"]])]
              } else {
                  multidata.info[[ p(ss) ]][["events.info"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- list(
                      event = ee,chrom=cc,
                      start = ddax$Start_bp[[ zzix.curr[[1]] ]],
                      end   = ddax$End_bp[[zzix.last]],
                      clonal = is.clonal,
                      idx    = p(zzix.curr,collapse=","),
                      N = length(zzix.curr)
                  )
                #multidata.info[[ p(ss) ]][["events.info"]][[length(multidata.info[[ p(ss) ]][["events.info"]])]] <- 
                #            p(multidata.info[[ p(ss) ]][["events.info"]][[length(multidata.info[[ p(ss) ]][["events.info"]])]],ddax$End_bp[[zzix.last]],":",clonality.string)
                  multidata.info[[ p(ss) ]][["events.idx"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- zzix.curr;
                  multidata.info[[ p(ss) ]][["events.data"]][[ length(multidata.info[[ p(ss) ]][["events.data"]])+1  ]] <- ddax[zzix.curr,,drop=F]
              }

              if(any( dda$event.type[idx[i:length(i)]] == ee )){
                if(DEBUG.MODE){ message("BREAK-step 2") }
                sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]+1;
                zzix <- which( ddax$rownum == dda$rownum[ idx[i] ])
                ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("START[",ee,"]"))
                if(DEBUG.MODE){ message("notes[i=",i,"]/[idx=",idx[i],"]/[ddax=",(zzix),"]=START") }
                multidata.info[[ p(ss) ]][["events"]] <- c(multidata.info[[ p(ss) ]][["events"]],p(ee,":chr",cc))
              }
            }  
          }
          zzix <- which( ddax$rownum == dda$rownum[ idx[length(idx)] ])
          ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("END[",ee,"]"))

              idx.prev   <- 1:length(idx);
              zzix.prev  <- which( ddax$rownum %in% dda$rownum[ idx[idx.prev] ] )
              notes.prev <- ddax$notes[ zzix.prev ];
              ixx.start  <- max( which(grepl(p("START[",ee,"]"),notes.prev,fixed=T) ))
              idx.curr   <- ixx.start:length(idx);
              zzix.curr  <- which( ddax$rownum %in% dda$rownum[ idx[idx.curr] ])
              zzix.last  <- zzix.curr[[length(zzix.curr)]]
              is.clonal  <- any( f.na( ddax$event.clonality[zzix.curr] == "clonal" ))
              is.event   <- any( f.na( ddax$event.type[zzix.curr] == ee ))
              clonality.string <- ifelse(is.clonal,"clonal","nonclonal")
              if(is.clonal){
                sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- sampCounts[[ss]][[p("mrgCT.clonal.",ee)]]+1;
                multidata.info[[ p(ss) ]][["events.clonal"]] <- c(multidata.info[[ p(ss) ]][["events.clonal"]],p(ee,":chr",cc))
              }
              if(is.event){
                  multidata.info[[ p(ss) ]][["events.info"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- list(
                      event = ee,chrom=cc,
                      start = ddax$Start_bp[[ zzix.curr[[1]] ]],
                      end   = ddax$End_bp[[zzix.last]],
                      clonal = is.clonal,
                      idx    = p(zzix.curr,collapse=","),
                      N = length(zzix.curr)
                  )
                #multidata.info[[ p(ss) ]][["events.info"]][[length(multidata.info[[ p(ss) ]][["events.info"]])]] <- 
                #            p(multidata.info[[ p(ss) ]][["events.info"]][[length(multidata.info[[ p(ss) ]][["events.info"]])]],ddax$End_bp[[zzix.last]],":",clonality.string)
                  multidata.info[[ p(ss) ]][["events.idx"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- zzix.curr;
                  multidata.info[[ p(ss) ]][["events.data"]][[ length(multidata.info[[ p(ss) ]][["events.data"]])+1  ]] <- ddax[zzix.curr,,drop=F]
              }
          if(DEBUG.MODE){ message("notes[i=",length(idx),"]/[idx=",idx[length(idx)],"]/[ddax=",(zzix),"]=END") }
          multidata.info[[ p(ss) ]][["data"]] <- ddax;
        }
      }
    }



  }
  sampCounts[[ss]][[p("mrgCT.ALL")]] <- 0;
  for(ee in eventlist ){
    sampCounts[[ss]][[p("mrgCT.ALL")]] <- sampCounts[[ss]][[p("mrgCT.ALL")]] + sampCounts[[ss]][[p("mrgCT.",ee)]]
  }

}

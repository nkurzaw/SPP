#' @import dplyr
annotateDataList <- function(dataList, geneNameVar, configLong,
                             intensityStr, fcStr, minQupm){
    # internal function to annotate list of 2D-TPP data subtables with
    # information from config table
    channel <- signal <- Temperature <- RefCol <- label <- 
        conc <- unique_ID <- spread_var <- NULL
    
    combinedTab <- bind_rows(lapply(dataList, function(dat){
        datLong <- dat %>% tbl_df() %>%
            gather(channel, signal, matches(intensityStr), 
                   matches(fcStr)) %>%
            mutate(label = gsub(fcStr, "", gsub(intensityStr, "", 
                                                channel))) %>%
            left_join(configLong %>% 
                          dplyr::select(Temperature, RefCol, label, 
                                        Experiment, conc, replicate),
                      by = c("temperature" = "Temperature", "label",
                             "experiment" = "Experiment")) %>% 
            group_by_(geneNameVar, "conc", "replicate") %>% 
            filter(qupm == max(qupm)) %>% 
            ungroup %>% 
            filter(qupm >= minQupm)
            
    }))
    return(combinedTab)
}

#' @import vsn
#' @importFrom tidyr spread
vsnNormalizeData <- function(datIn, idVar, intensityStr){
    datSpread <- datIn %>% 
        mutate(channel_rep = paste(channel, replicate, sep = "__")) %>% 
        dplyr::select(representative, clustername, channel_rep, signal) %>% 
        spread(channel_rep, signal)
    
    vsn_fit <- vsn::vsn2(as.matrix(datSpread[,-c(1,2)]))
    vsnNorm <- vsn::predict(vsn_fit, as.matrix(datSpread[,-c(1,2)]))
    datNorm <- cbind(datSpread[,c(1,2)], 2^vsnNorm)
    datNormLong <- as_tibble(datNorm) %>% 
        gather(channel_rep, norm_signal, -representative, -clustername) %>% 
        separate(channel_rep, c("channel", "replicate"), sep = "__") %>% 
        mutate(replicate = as.factor(replicate))
    
    out_df <- left_join(datIn, datNormLong, 
              by = c("representative", "clustername", 
                     "channel", "replicate")) %>% 
        mutate(log2_signal = log2(norm_signal))
    return(out_df)
}

#' Import SPP dataset using a config table
#' 
#' @param configTable character string of a file path to a config table
#' @param data possible list of datasets from different MS runs 
#' corresponding to a 2D-TPP dataset, circumvents loading datasets 
#' referencend in config table, default is NULL
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#' @param intensityStr character string indicating which columns contain 
#'   raw intensities measurements
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the suffix \code{fcStr} 
#'   will be regarded as containing fold change values.
#' @param naStrs character vector indicating missing values in the data table. 
#'   When reading data from file, this value will be passed on to the argument 
#'   \code{na.strings} in function \code{read.delim}.
#' @param qualColName character string indicating which column can be used for 
#'   additional quality criteria when deciding between different non-unique 
#'   protein identifiers.
#' @param medianNormalizeFC perform median normalization (default: TRUE).
#' @param addCol character string indicating additional column to import
#' @param filterContaminants boolean variable indicating whether data 
#' should be filtered to exclude contaminants (default: TRUE).
#' @param nonZeroCols column like default qssm that should be imported and
#' requested to be non-zero in analyzed data
#' @param geneNameVar character string of the column name that describes
#' the gene name of a given protein in the raw data files
#' @param concFactor numeric value that indicates how concentrations need to 
#' be adjusted to yield total unit e.g. default mmol - 1e6
#' 
#' @return tidy data frame representing a 2D-TPP dataset
#' 
#' @examples 
#' data("config_tab")
#' data("raw_dat_list")
#' import_df <- import2dDataset(configTable = config_tab, 
#'                              data = raw_dat_list,
#'                              idVar = "protein_id",
#'                              intensityStr = "signal_sum_",
#'                              fcStr = "rel_fc_",
#'                              nonZeroCols = "qusm",
#'                              geneNameVar = "gene_name",
#'                              addCol = NULL,
#'                              qualColName = "qupm",
#'                              naStrs = c("NA", "n/d", "NaN"),
#'                              concFactor = 1e6,
#'                              medianNormalizeFC = TRUE,
#'                              filterContaminants = TRUE)
#' 
#' @export
#' 
#' @import TPP2D
#' @import dplyr
#' @import tidyr
#' @import vsn
importSppDataset <- function(configTable, data,
                            idVar = "representative",
                            intensityStr = "sumionarea_protein_",
                            nonZeroCols = "qssm",
                            geneNameVar = "clustername",
                            addCol = NULL,
                            qualColName = "qupm",
                            naStrs = c("NA", "n/d", "NaN"),
                            concFactor = 1e6,
                            minQupm = 2,
                            vsnNormalize = TRUE,
                            filterContaminants = TRUE){
    
    configWide <- TPP2D:::TPP_importCheckConfigTable(
        infoTable = configTable, type = "2D")
    configLong <- TPP2D:::configWide2Long(configWide = configWide) %>% 
        mutate(replicate = as.factor(
            factor(dense_rank(Experiment), 
                   levels = seq_len(length(unique(Experiment))))))
    
    dataList <- TPP2D:::import2dMain(configTable = configWide,
                             data = data,
                             idVar = idVar,
                             fcStr = NULL,
                             addCol = c(geneNameVar, addCol),
                             naStrs = naStrs,
                             intensityStr = intensityStr,
                             nonZeroCols = nonZeroCols,
                             qualColName = qualColName)
    
    dataLong <- annotateDataList(dataList = dataList,
                                 geneNameVar = geneNameVar,
                                 configLong = configLong,
                                 intensityStr = intensityStr,
                                 fcStr = fcStr,
                                 minQupm = minQupm)
    
    dataRenamed <- TPP2D:::renameColumns(dataLong = dataLong,
                                         idVar = idVar,
                                         geneNameVar = geneNameVar)
    
    if(filterContaminants){
        dataRenamed <- TPP2D:::filterOutContaminants(dataRenamed)
    }
    if(vsnNormalize){
        dataNorm <- vsnNormalizeData(dataRenamed, idVar = idVar, 
                                     intensityStr = intensityStr)
    }else{
        dataNorm <- dataRenamed %>% 
            mutate(log2_signal = log2(signal))
    }
    
    return(dataNorm)
}

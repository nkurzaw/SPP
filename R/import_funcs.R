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
#' @import vsn
importSppDataset <- function(configTable, data,
                            idVar = "representative",
                            intensityStr = "sumionarea_protein_",
                            fcStr = "rel_fc_protein_",
                            nonZeroCols = "qssm",
                            geneNameVar = "clustername",
                            addCol = NULL,
                            qualColName = "qupm",
                            naStrs = c("NA", "n/d", "NaN"),
                            concFactor = 1e6,
                            vsnNormalize = TRUE,
                            filterContaminants = TRUE){
    
    configWide <- TPP2D:::TPP_importCheckConfigTable(
        infoTable = configTable, type = "2D")
    configLong <- TPP2D:::configWide2Long(configWide = configWide) %>% 
        mutate(replicate = as.factor(
            factor(dense_rank(Experiment), 
                   levels = seq_len(length(unique(Experiment))))))
    
    dataList <- import2dMain(configTable = configWide,
                             data = data,
                             idVar = idVar,
                             fcStr = fcStr,
                             addCol = c(geneNameVar, addCol),
                             naStrs = naStrs,
                             intensityStr = intensityStr,
                             nonZeroCols = nonZeroCols,
                             qualColName = qualColName)
    
    dataLong <- annotateDataList(dataList = dataList,
                                 geneNameVar = geneNameVar,
                                 configLong = configLong,
                                 intensityStr = intensityStr,
                                 fcStr = fcStr)
    
    dataRatioChecked <- checkRatioRef(dataLong, idVar = idVar,
                                      concFactor = concFactor)
    
    if(medianNormalizeFC){
        message("Median normalizing fold changes...")
        dataNorm <- medianNormalizeRatios(dataRatioChecked)
    }else{
        dataNorm <- dataRatioChecked
    }
    
    dataOut <- renameColumns(dataNorm,
                             idVar = idVar,
                             geneNameVar = geneNameVar)
    
    if(filterContaminants){
        dataOut <- filterOutContaminants(dataOut)
    }
    
    return(dataOut)
}
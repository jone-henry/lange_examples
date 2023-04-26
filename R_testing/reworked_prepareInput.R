prepareInput <- function(
    data_file_w_path,
    meta_file_w_path,
    metaColName
){
  
  # Sanity check - do the files exist/are the paths correct
  # If they are not, then it will print a descriptive error message
  data_check = file.exists(data_file_w_path)
  meta_check = file.exists(meta_file_w_path)
  if(data_check == FALSE){
    stop("! Data file does not exist and/or path to file in incorrect ! - verify and re-try")}
  if(meta_check == FALSE){
    stop("! Metadata file does not exist and/or path to file in incorrect ! - verify and re-try")}
  
  # Sanity check - is there actually anything in the files
  # If they are empty it will print a descriptive error message
  # Note size is given in bytes
  data_size_check = file.info(data_file_w_path)$size
  meta_size_check = file.info(meta_file_w_path)$size
  if(data_size_check == 0){
    stop("! Data file is empty ! - verify file contents and re-try")}
  if(meta_size_check == 0){
    stop("! Metadata file is 0 bytes ! - verify file contents and re-try")}
  
  # Open Datasets
  data_df     = data.frame(read.csv(data_file_w_path, header=T, check.names=F, sep=',', 
                                    na.strings=c('',' ','NA','NaN','na','nan','NAN','Nan','Filtered')))
  metadata_df = data.frame(read.csv(meta_file_w_path, header=T, check.names=F, sep=',', 
                                    na.strings=c('',' ','NA','NaN','na','nan','NAN','Nan')))
  # Quick visual check of the data
  # Headers +2 lines of data will be printed
  cat("Input data headers + first two lines of data:\n")
  print(head(data_df,2))
  cat("Dataframe dimensions:", dim(data_df))
  cat('\n')
                           
  # Now that the data is in, some more sanity checks
  # Check that there is actually data. Note this is different from the above
  # An errant 'cat' command or rogue redirect can make an empty file
  # But if there is an error in creating the file elsewhere, then it may make
  #   a file with just headers and nothing else in it
  # Checking for <1 rows of data
  if(nrow(data_df) < 1){
    stop("! Data file does not appear to contain any data ! - verify file contents and re-try")}
  if(nrow(metadata_df) < 1){
    stop("! Metadata file does not appear to contain any data ! - verify file contents and re-try")}
  
  # Check and see if the provided metaColName is actually in the file
  if(metaColName %in% colnames(metadata_df) == FALSE)
    stop("! Metadata file does NOT contain the column title passed in ! - verify input data and re-try")
  
  # Check for NA in metadata column names
  # If any are found, inform the user and remove
  if(any(is.na(metadata_df[, metaColName])) == TRUE){
    print("Note: NA rows have been found in the metadata file and will be removed", quote=FALSE)
    metadata_df = data.frame(metadata_df[!is.na(metadata_df[, metaColName]), ])
    print('Data import WILL continue', quote=FALSE)}
  
  col2select <- intersect(metadata_df[, metaColName], colnames(data_df))
  # Warn the user when the number of columns selected from data_df is less
  #   than the number of anticipated columns according to meta_df[, metaColName]
  if (length(col2select) < length(metadata_df[, metaColName])){
    print("Note 1: Not all columns indicated in the metadata file could be found in the data file", quote=FALSE)
    print("Note 2: Data import WILL continue", quote=FALSE)
    print("If this is unexpected, verify input data and retry", quote=FALSE)}
  # Also warn the user if no columns intersect
  intersection_check = length(col2select)
  if(intersection_check==0){
    stop("! No intersecting columns were found between the metadata file and the data file ! - verify input data and re-try")}
  
  # Check for duplicate proteins in the first column of the data file
  # Warn user of duplicates and remove doubles, triples, etc.
  if(any(duplicated(data_df[,1])==TRUE)){
    row_names        = data_df[,1]
    duplicate_rnames = row_names[duplicated(row_names)]
    print("Note: The following protein(s) were found to have 2 (or more) times", quote=FALSE)
    print(unique(duplicate_rnames))
    print('Only the first occurence has been kept', quote=FALSE)
    print('If this was unexpected, please review your data and retry', quote=FALSE)
    data_df = data_df[!duplicated(data_df[,1]), ]}
  
  # rename rows to index values 
  rownames(data_df) <- 1:nrow(data_df)

  # Convert elements to numeric in the quantitative part of the data
  quant_data <- sapply(data_df[, col2select], as.numeric)
  # Add rownames to quant.data
  rownames(quant_data) <- rownames(data_df)
  
  # Simultaneously remove both columns and rows where quantitative data is COMPLETELY missing
  quant_data_processed = data.frame(quant_data[!rowSums(is.na(quant_data)) == ncol(quant_data),
                                               !colSums(is.na(quant_data)) == nrow(quant_data)])
  if(identical(dim(quant_data), dim(quant_data_processed)) == FALSE){
    print('Note: Some columns and/or rows have been removed', quote=FALSE)
    cat('Old dimensions: ', dim(quant_data))
    cat('\n')
    cat('New dimensions: ', dim(quant_data_processed))
    cat('\n')
    cat('Row indices removed (if applicable):', names(which(rowSums(is.na(quant_data)) == ncol(quant_data))))
    cat('\n')
    cat('Column(s) removed (if applicable):',   names(which(colSums(is.na(quant_data)) == nrow(quant_data))))
    cat('\n')}
    
  # Make sure the rownames are consistent with annotation and quantitative data
  data_df = data_df[rownames(quant_data_processed), colnames(quant_data_processed)]
  # Note that this will return the data but the numeric fields may not all be actually numeric (since we did apply as.numeric only to quant.data for the subset columns)
  cat("Data successfully imported and processed")
  return(data_df)
}

# Reading in known good files
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example.csv',
                   'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                   'Columns')

# Checking how it handles missing files
# data
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\000.csv',
                   'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                   'Columns')
# metadata
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\000.csv',
                    'Columns')

# Checking how it handles empty files
# data
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example_empty.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                    'Columns')
# metadata
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2_empty.csv',
                    'Columns')

# Checking how it handles files containing only headers
# data
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example_headers_only.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                    'Columns')
# metadata
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2_headers_only.csv',
                    'Columns')

# Check how it handles having an incorrect 3rd variable
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                    'ABC')

# Check how it handles when the data file is missing one of the columns indicated in the metadata file
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example_missing_column.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                    'Columns')

# Check how it handles when the data file has NONE of the columns indicated in the metadata file
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example_bad_columns.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                    'Columns')
# opposite of above - bad metadata column names
# Check how it handles when the metadata file gives incorrect column names
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2_bad_columns.csv',
                    'Columns')

# Check how it handles when the metadata file has a blank in one of the rows in the "Columns" column
# Should read in as NA and inform the user
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2_NA.csv',
                    'Columns')

# Check how it handles when the data file has a blank in one of the cells
# Should read in as NA and proceed as normal (no warning)
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example_blank_cell.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                    'Columns')

# Check how it handles a duplicate protein
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example_duplicated_row.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                    'Columns')

# Check how it handles a triplicate protein
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example_triplicate_row.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                    'Columns')

# Check how it handles a row of NA
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example_one_row_NA.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                    'Columns')

# Check how it handles a col of NA
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example_one_col_NA.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                    'Columns')

# Check how it handles one col + 1 row of NA
test = prepareInput('C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\bmif-Example_one_col_one_row_NA.csv',
                    'C:\\Users\\John\\PycharmProjects\\Lange\\R_testing\\Metadata-Example-2.csv',
                    'Columns')

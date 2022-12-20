library(h2o)

#Function to removing everything but the train and test dataset from the h2o memory
rm_h2o_elements <- function(train_test_keys){
    #' Function to remove all elements in H2O cluster but the training and testing ones
    #'
    #' @param train_test_keys     list of h2o keys that just refer to the training and testing datasets, i.e. things that have to be kept
    #'
    #' @returns                   does not return anything but cleans the h2o cluster memory

  key_list <- h2o::h2o.ls()$key
  #Get list of keys to be deleted
  delete_list <- unlist(key_list[!(key_list %in% train_test_keys) & grepl("RTMP", key_list, fixed = TRUE) == FALSE])
  #Remove those elements
  h2o::h2o.rm(delete_list)
}
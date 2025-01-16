# Load required library
library(dplyr)
library(tidyverse)
library(stringr)

#Load the QC limit file
limit_df <- read.csv("~/Downloads/Limit_QC.csv", header = TRUE, check.names = FALSE)

# Define function to extract lines after finding specific "DataType:" conditions
extract_sections <- function(file_path, vaccine_input) {
  # Read the CSV file, skipping the first 51 lines
  data <- read.csv(file_path, header = FALSE, stringsAsFactors = FALSE, fill = TRUE, skip = 51)
  
  # Helper function to extract rows between a start and end criteria
  extract_data <- function(start_text) {
    # Find the row that contains the specified "DataType:" text
    start_row <- which(data$V1 == "DataType:" & data$V2 == start_text)
    
    # Check if the start criteria is found
    if (length(start_row) == 0) {
      message(paste("The line with 'DataType:' and", start_text, "was not found in the file."))
      return(NULL)
    }
    
    # Find the end row where a line starts with whitespace or another "DataType:"
    end_row <- start_row + 1
    while (end_row <= nrow(data) && 
           !grepl("^\\s*$", data[end_row, 1]) && 
           !(data[end_row, 1] == "DataType:")) {
      end_row <- end_row + 1
    }
    
    # Extract rows between start and end criteria
    extracted_data <- data[(start_row + 1):(end_row - 1), ]
    
    # Use first row as headers for columns starting with "V" and remove that row
    colnames(extracted_data) <- ifelse(grepl("^V", colnames(extracted_data)), extracted_data[1, ], colnames(extracted_data))
    extracted_data <- extracted_data[-1, ]  # Remove the first row
    return(extracted_data)
  }
  
  # Extract data for each section
  result_data <- extract_data("Result")
  r2_data <- extract_data("R^2")
  count_data <- extract_data("Count")
  
  # Add Group, Vaccine, Sample_ID, and Stimulus columns to result_data, and filter out "Background" or "Standard" samples
  if (!is.null(result_data)) {
    result_data <- result_data %>%
      mutate(
        Group = str_extract(Sample, "^[0-9]+"),                   # Extract the first number for Group
        Vaccine = vaccine_input[Group],                           # Map Vaccine using user input for each Group
        Sample_ID = str_extract(Sample, "^[0-9]+[A-Za-z]"),           # Extract the Sample ID
        Stimulus = "Placeholder"
      ) %>%
      filter(!grepl("Background|Standard|STD|BLANK", Sample, ignore.case = TRUE)) %>%  # Exclude rows with "Background" or "Standard"
      select(1:2, Group, Vaccine, Sample_ID, Stimulus, everything())  # Move new columns to 3rd-6th position
  }
  count_data <- count_data %>% filter(!grepl("Background|Standard|STD|BLANK", Sample, ignore.case = TRUE))
  # Store each result in a named list
  list(
    Result = result_data,
    R2 = r2_data,
    Count = count_data
  )
}


# Usage example
file_path <- "~/Downloads/YpKup_LTA1_TCSups_24Oct24replay.csv"  # Replace with the path to your CSV file



# User input for vaccines by group
vaccine_input <- list("A" = "Vaccine1", "B" = "Vaccine2", "C" = "Vaccine3")  # Update with actual group-vaccine mappings

# Extract data frames for "Result", "R^2", and "Count"
data_frames <- extract_sections(file_path, vaccine_input)

# Access each data frame
result_df <- data_frames$Result
r2_df <- data_frames$R2
count_df <- data_frames$Count

############ remove total events column and other unucessary tail columns
result_df <- result_df[1:(length(result_df)-1)]
count_df <- count_df[1:(length(count_df)-1)]
r2_df <- r2_df[1:(length(r2_df)-2)]


######## Find unique groups
unique_groups <- unique(result_df$Group)

# Initialize an empty list to store user inputs for each group
vaccine_inputs <- list()

# Loop through each unique group and ask user for input
for (group in unique_groups) {
  prompt <- paste("Enter vaccine type for Group", group, ": ")
  
  # Get user input and store it in the list with the group as the key
  vaccine_inputs[[group]] <- readline(prompt)
}

# Apply the user inputs to the Vaccine column based on the Group
result_df$Vaccine <- sapply(result_df$Group, function(g) vaccine_inputs[[g]])






# Extract the unique last numbers directly from the Sample column
unique_groups_stim <- unique(str_extract(result_df$Sample, "[0-9]+$"))

# Initialize an empty list to store user inputs for each unique last number
stimulus_inputs <- list()

# Loop through each unique group and ask user for input
for (stimgroup in unique_groups_stim) {
  prompt <- paste("Input stimulus for Sample Group X", stimgroup, ": ", sep = "")
  
  # Get user input and store it in the list with the group as the key
  stimulus_inputs[[stimgroup]] <- readline(prompt)
}

# Apply the user inputs to the Stimulus column based on the last character of each sample
result_df$Stimulus <- sapply(result_df$Sample, function(sample) {
  stimgroup <- substr(sample, nchar(sample), nchar(sample))
  stimulus_inputs[[stimgroup]]
})

#Logic to remove any blanks from stimulus (eg. any empty wells/Samples) 
# could be more elegant on the front end and ask user for total Samples used
result_df <- result_df[!result_df$Stimulus == "",]



# Find the starting index for the "Stimulus" column and identify all columns that follow it
start_col <- which(names(result_df) == "Stimulus") + 1  # +1 to start after "Stimulus"
target_columns <- names(result_df)[start_col:ncol(result_df)]  # Extract names of columns from "G-CSF" onward

## Clean columns from "G-CSF" onward: remove "<" or ">" symbols and convert to numeric
result_df[target_columns] <- lapply(result_df[target_columns], function(column) {
  as.numeric(gsub("[<>]", "", column))
})

#Ensure target columns exist in both count_df and result_df
if(all(target_columns %in% names(count_df))) {
  
  # Loop over each target column to check for values < 35 in count_df
  for (col in target_columns) {
    # Identify rows in count_df where the value is < 35 for the current column
    rows_to_replace <- which(as.numeric(count_df[[col]]) < 35)
    
    # Replace corresponding values in result_df with NA
    result_df[rows_to_replace, col] <- NA
  }
  
  # Display modified result_df
  print("Modified 'result_df' with values replaced by NA based on count_df criteria:")
  #print(result_df)
} else {
  warning("Target columns from count_df are not all present in result_df.")
}




# Extract the ULOQ and LLOQ rows based on their labels in the first column
uloq_row <- which(limit_df[, 1] == "ULOQ")
lloq_row <- which(limit_df[, 1] == "LLOQ")

# Convert the column names to match those in result_df
limit_columns <- names(limit_df)[-1]  # Exclude the first column, which contains "ULOQ" and "LLOQ" labels

# Loop through each target column in result_df
for (col in target_columns) {
  # Check if the column in result_df has a corresponding ULOQ and LLOQ column in limits
  if (col %in% limit_columns) {
    # Find the index of the column in the limits data to match the correct ULOQ and LLOQ values
    limit_index <- which(limit_columns == col)
    
    # Retrieve ULOQ and LLOQ values for this specific analyte
    uloq <- as.numeric(limit_df[uloq_row, limit_index + 1])  # Adjust for the excluded first column
    lloq <- as.numeric(limit_df[lloq_row, limit_index + 1])  # Adjust for the excluded first column
    
    # Replace values in result_df that exceed ULOQ with the ULOQ value
    result_df[[col]][result_df[[col]] > uloq] <- uloq
    
    # Replace values in result_df that are below LLOQ with the LLOQ value
    result_df[[col]][result_df[[col]] < lloq] <- lloq
  }
}

# Define a function to calculate the geometric mean
geo_mean <- function(x) {
  exp(mean(log(x[x > 0]), na.rm = TRUE))  # Handle zeros by excluding them from the log transformation
}

# Aggregate technical replicates (rows with identical `Sample` values) by calculating the mean
result_df <- result_df %>%
  group_by(Sample) %>%
  summarize(across(all_of(target_columns), ~ mean(as.numeric(.), na.rm = TRUE)), .groups = "drop") %>%
  left_join(result_df %>% select(Sample, Group, Vaccine, Stimulus, Sample_ID) %>% distinct(), by = "Sample")

# Calculate geometric mean for each group and stimulus across columns for all target_columns (cytokines)
geometric_means_df <- result_df %>%
  group_by(Group, Stimulus) %>%
  summarize(
    Vaccine = first(Vaccine),
    Sample_ID = first(Sample_ID),
    across(all_of(target_columns), geo_mean), .groups = "drop")

# Optional: Rename columns back to the original names for easier comparison
#names(geometric_means_df)[3:ncol(geometric_means_df)] <- names(data)[start_col:ncol(data)]


####Should add could here to intake upper and lower limit QC file to alter values above/below those limits####


# Display the results
#print("Data for 'DataType: Result':")
#print(result_df)

#print("Data for 'DataType: R^2':")
#print(r2_df)

#print("Data for 'DataType: Count':")
#print(count_df)

write.csv(geometric_means_df, "~/Downloads/final_out.csv", row.names = FALSE)

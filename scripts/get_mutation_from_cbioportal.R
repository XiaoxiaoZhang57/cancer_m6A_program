# Complete Code for Extracting Genomic Positions of Cancer Mutation Sites
# This code retrieves mutation data from all cancer studies in cBioPortal database and extracts genomic position information

# 1. Install and load necessary packages --------------------------------------------------------
if (!require("cbioportalR")) {
    install.packages("cbioportalR")
}
if (!require("dplyr")) {
    install.packages("dplyr")
}
if (!require("purrr")) {
    install.packages("purrr")
}

library(cbioportalR)
library(dplyr)
library(purrr)

# 2. Set up cBioPortal database connection -------------------------------------------------
set_cbioportal_db("https://www.cbioportal.org/api")
print("cBioPortal database connection established")

# 3. Get all studies -------------------------------------------------------------
get_all_studies_info <- function() {
    print("Getting list of all studies...")
    all_studies <- available_studies()
    print(paste("Found", nrow(all_studies), "studies"))
    
    # Display study information
    print("Study data structure:")
    print(str(all_studies))
    
    print("First 10 studies:")
    print(all_studies[1:10, c("studyId", "name", "cancerTypeId")])
    
    return(all_studies)
}

# 4. Function to get mutation data -------------------------------------------------------
get_mutations_simple <- function(study_id) {
    tryCatch({
        cat("Processing study:", study_id, "\n")
        
        # Directly use study_id to get mutation data
        mutations <- get_mutations_by_study(study_id = study_id)
        
        if (is.null(mutations) || nrow(mutations) == 0) {
            cat("  Study", study_id, "has no mutation data\n")
            return(NULL)
        }
        
        # Add study information
        study_info <- all_studies %>% filter(studyId == study_id)
        mutations$study_id <- study_id
        mutations$cancer_type <- ifelse(nrow(study_info) > 0, study_info$cancerTypeId, "Unknown")
        mutations$study_name <- ifelse(nrow(study_info) > 0, study_info$name, "Unknown")
        
        cat("  Successfully retrieved", nrow(mutations), "mutations\n")
        return(mutations)
        
    }, error = function(e) {
        cat("  Error processing study", study_id, ":", e$message, "\n")
        return(NULL)
    })
}

# 5. Function to batch process all studies ---------------------------------------------------
process_all_studies <- function(studies_df, batch_size = 5, max_studies = NULL) {
    # If max_studies is specified, limit the number of studies to process
    if (!is.null(max_studies)) {
        studies_to_process <- studies_df %>% head(max_studies)
    } else {
        studies_to_process <- studies_df
    }
    
    total_studies <- nrow(studies_to_process)
    print(paste("Starting to process", total_studies, "studies, batch size:", batch_size))
    
    all_results <- list()
    processed_count <- 0
    error_count <- 0
    
    # Process studies in batches
    for (i in seq(1, total_studies, batch_size)) {
        batch_end <- min(i + batch_size - 1, total_studies)
        current_batch <- studies_to_process$studyId[i:batch_end]
        
        print(paste("Processing batch:", i, "to", batch_end, "/", total_studies))
        
        # Process current batch of studies
        batch_results <- map(current_batch, function(study_id) {
            result <- get_mutations_simple(study_id)
            processed_count <<- processed_count + 1
            return(result)
        })
        
        # Filter out NULL results
        valid_results <- batch_results %>% compact()
        
        if (length(valid_results) > 0) {
            # Combine results from current batch
            batch_df <- bind_rows(valid_results)
            all_results[[length(all_results) + 1]] <- batch_df
            
            print(paste("Batch completed, retrieved", nrow(batch_df), "mutations"))
        } else {
            print("No mutation data retrieved in current batch")
        }
        
        # Pause to avoid too frequent requests
        Sys.sleep(2)
    }
    
    # Combine results from all batches
    if (length(all_results) > 0) {
        final_df <- bind_rows(all_results)
        print(paste("All batches processed! Total of", nrow(final_df), "mutation records retrieved"))
        print(paste("Successfully processed", processed_count, "studies"))
        return(final_df)
    } else {
        print("No mutation data retrieved")
        return(NULL)
    }
}

# 6. Function to extract and save genomic position information -------------------------------------------
extract_and_save_genomic_positions <- function(mutations_df, output_prefix = "cancer_mutations") {
    if (is.null(mutations_df) || nrow(mutations_df) == 0) {
        print("No mutation data to process")
        return(NULL)
    }
    
    print("Starting to extract genomic position information...")
    
    # Define required columns
    required_cols <- c("chr", "startPosition", "endPosition", "referenceAllele", 
                      "variantAllele", "mutationType", "proteinChange", 
                      "hugoGeneSymbol", "entrezGeneId")
    
    # Check which columns exist
    available_cols <- required_cols[required_cols %in% names(mutations_df)]
    print(paste("Available genomic position columns:", paste(available_cols, collapse = ", ")))
    
    # Extract genomic position information
    genomic_data <- mutations_df %>%
        select(all_of(c(available_cols, "study_id", "cancer_type", "study_name"))) %>%
        distinct() %>%
        filter(!is.na(chr) & !is.na(startPosition))
    
    # Ensure correct data types
    genomic_data$chr <- as.character(genomic_data$chr)
    genomic_data$startPosition <- as.numeric(genomic_data$startPosition)
    genomic_data$endPosition <- as.numeric(genomic_data$endPosition)
    genomic_data$hugoGeneSymbol <- as.character(genomic_data$hugoGeneSymbol)
    
    print(paste("Extracted", nrow(genomic_data), "unique genomic positions"))
    
    # Save complete mutation data
    complete_output_file <- paste0(output_prefix, "_complete.csv")
    write.csv(mutations_df, complete_output_file, row.names = FALSE)
    print(paste("Complete mutation data saved to:", complete_output_file))
    
    # Save genomic position data
    genomic_output_file <- paste0(output_prefix, "_genomic_positions.csv")
    write.csv(genomic_data, genomic_output_file, row.names = FALSE)
    print(paste("Genomic position data saved to:", genomic_output_file))
    
    # Create BED format file (for genome browser)
    bed_output_file <- paste0(output_prefix, "_genome_browser.bed")
    bed_data <- genomic_data %>%
        mutate(
            chrom = ifelse(grepl("^chr", chr), chr, paste0("chr", chr)),
            name = paste(hugoGeneSymbol, 
                        ifelse(!is.na(proteinChange) & proteinChange != "", proteinChange, "NA"), 
                        cancer_type, study_id, sep = "|"),
            score = 1000,
            strand = "."
        ) %>%
        select(chrom, startPosition, endPosition, name, score, strand)
    
    write.table(bed_data, bed_output_file, 
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    print(paste("BED format file saved to:", bed_output_file))
    
    # Create summary file grouped by cancer type
    cancer_summary <- genomic_data %>%
        group_by(cancer_type) %>%
        summarise(
            mutation_count = n(),
            unique_genes = n_distinct(hugoGeneSymbol),
            studies = n_distinct(study_id),
            .groups = 'drop'
        ) %>%
        arrange(desc(mutation_count))
    
    summary_output_file <- paste0(output_prefix, "_cancer_summary.csv")
    write.csv(cancer_summary, summary_output_file, row.names = FALSE)
    print(paste("Cancer type summary saved to:", summary_output_file))
    
    # Create gene summary file
    gene_summary <- genomic_data %>%
        group_by(hugoGeneSymbol, chr) %>%
        summarise(
            mutation_count = n(),
            cancer_types = n_distinct(cancer_type),
            studies = n_distinct(study_id),
            min_position = min(startPosition, na.rm = TRUE),
            max_position = max(endPosition, na.rm = TRUE),
            .groups = 'drop'
        ) %>%
        arrange(desc(mutation_count))
    
    gene_output_file <- paste0(output_prefix, "_gene_summary.csv")
    write.csv(gene_summary, gene_output_file, row.names = FALSE)
    print(paste("Gene summary saved to:", gene_output_file))
    
    return(genomic_data)
}

# 7. Function to save separate files by cancer type ---------------------------------------------
save_by_cancer_type <- function(genomic_data, output_dir = "cancer_type_files") {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }
    
    cancer_types <- unique(genomic_data$cancer_type)
    print(paste("Creating separate files for", length(cancer_types), "cancer types..."))
    
    for (cancer in cancer_types) {
        cancer_data <- genomic_data %>% filter(cancer_type == cancer)
        if (nrow(cancer_data) > 0) {
            # Create safe filename
            safe_cancer_name <- gsub("[^a-zA-Z0-9]", "_", cancer)
            filename <- file.path(output_dir, paste0("mutations_", safe_cancer_name, ".csv"))
            write.csv(cancer_data, filename, row.names = FALSE)
        }
    }
    
    print(paste("Mutation data for each cancer type saved to", output_dir, "directory"))
}

# 8. Function to generate statistical report -------------------------------------------------------
generate_statistics_report <- function(genomic_data, mutations_df) {
    print("Generating statistical report...")
    
    # Ensure correct data types
    genomic_data$chr <- as.character(genomic_data$chr)
    genomic_data$cancer_type <- as.character(genomic_data$cancer_type)
    genomic_data$hugoGeneSymbol <- as.character(genomic_data$hugoGeneSymbol)
    if ("mutationType" %in% names(genomic_data)) {
        genomic_data$mutationType <- as.character(genomic_data$mutationType)
    }
    
    # Basic statistics
    total_mutations <- nrow(mutations_df)
    unique_genomic_positions <- nrow(genomic_data)
    unique_genes <- length(unique(genomic_data$hugoGeneSymbol))
    unique_cancer_types <- length(unique(genomic_data$cancer_type))
    unique_studies <- length(unique(genomic_data$study_id))
    
    # Use table function instead of count to avoid data type issues
    chr_distribution <- as.data.frame(table(genomic_data$chr), stringsAsFactors = FALSE)
    colnames(chr_distribution) <- c("chr", "n")
    chr_distribution <- chr_distribution %>% arrange(desc(n))
    
    # Cancer type distribution
    cancer_distribution <- as.data.frame(table(genomic_data$cancer_type), stringsAsFactors = FALSE)
    colnames(cancer_distribution) <- c("cancer_type", "n")
    cancer_distribution <- cancer_distribution %>% arrange(desc(n))
    
    # Gene distribution
    gene_distribution <- as.data.frame(table(genomic_data$hugoGeneSymbol), stringsAsFactors = FALSE)
    colnames(gene_distribution) <- c("hugoGeneSymbol", "n")
    gene_distribution <- gene_distribution %>% arrange(desc(n))
    
    # Mutation type distribution (if exists)
    if ("mutationType" %in% names(genomic_data)) {
        mutation_type_distribution <- as.data.frame(table(genomic_data$mutationType), stringsAsFactors = FALSE)
        colnames(mutation_type_distribution) <- c("mutationType", "n")
        mutation_type_distribution <- mutation_type_distribution %>% arrange(desc(n))
    } else {
        mutation_type_distribution <- NULL
    }
    
    # Create report
    report <- list(
        summary = list(
            total_mutations = total_mutations,
            unique_genomic_positions = unique_genomic_positions,
            unique_genes = unique_genes,
            unique_cancer_types = unique_cancer_types,
            unique_studies = unique_studies
        ),
        chromosome_distribution = chr_distribution,
        cancer_type_distribution = cancer_distribution,
        gene_distribution = gene_distribution,
        mutation_type_distribution = mutation_type_distribution
    )
    
    # Save report as CSV
    summary_df <- data.frame(
        Metric = names(report$summary),
        Value = unlist(report$summary)
    )
    write.csv(summary_df, "mutation_statistics_summary.csv", row.names = FALSE)
    
    # Save distribution data
    write.csv(chr_distribution, "chromosome_distribution.csv", row.names = FALSE)
    write.csv(cancer_distribution, "cancer_type_distribution.csv", row.names = FALSE)
    write.csv(gene_distribution, "gene_distribution.csv", row.names = FALSE)
    
    if (!is.null(mutation_type_distribution)) {
        write.csv(mutation_type_distribution, "mutation_type_distribution.csv", row.names = FALSE)
    }
    
    print("Statistical report generated and saved")
    return(report)
}

# 9. Main execution function ---------------------------------------------------------------
main <- function(max_studies = NULL, batch_size = 5) {
    print("=== Starting batch processing of all cancer studies ===")
    
    # 1. Get all studies
    print("Step 1: Getting study list")
    all_studies <- get_all_studies_info()
    
    # 2. Batch process all studies
    print("Step 2: Batch retrieving mutation data")
    all_mutations <- process_all_studies(all_studies, 
                                        batch_size = batch_size, 
                                        max_studies = max_studies)
    
    if (is.null(all_mutations)) {
        print("No mutation data retrieved, program ending")
        return(NULL)
    }
    
    # 3. Extract and save genomic position information
    print("Step 3: Extracting and saving genomic position information")
    genomic_positions <- extract_and_save_genomic_positions(all_mutations, "all_cancer_mutations")
    
    if (is.null(genomic_positions)) {
        print("No genomic position data to save")
        return(NULL)
    }
    
    # 4. Save separate files by cancer type
    print("Step 4: Saving files by cancer type")
    save_by_cancer_type(genomic_positions)
    
    # 5. Generate statistical report
    print("Step 5: Generating statistical report")
    stats_report <- generate_statistics_report(genomic_positions, all_mutations)
    
    # 6. Display final summary
    print("=== Processing completed ===")
    print(paste("Total mutation records:", stats_report$summary$total_mutations))
    print(paste("Unique genomic positions:", stats_report$summary$unique_genomic_positions))
    print(paste("Number of genes involved:", stats_report$summary$unique_genes))
    print(paste("Number of cancer types:", stats_report$summary$unique_cancer_types))
    print(paste("Number of studies:", stats_report$summary$unique_studies))
    
    # 7. Generate final summary report
    generate_final_summary(stats_report, genomic_positions)
    
    return(list(
        mutations = all_mutations,
        genomic_positions = genomic_positions,
        statistics = stats_report
    ))
}

# 10. Function to generate final summary report --------------------------------------------------
generate_final_summary <- function(stats_report, genomic_positions) {
    print("Generating final summary report...")
    
    summary_text <- c(
        "Final Summary Report of Cancer Mutation Database",
        paste("Generated at:", Sys.time()),
        paste("Data source: cBioPortal database (https://www.cbioportal.org/api)"),
        "",
        "=== Data Scale Summary ===",
        paste("Total mutation records:", format(stats_report$summary$total_mutations, big.mark = ",")),
        paste("Unique genomic positions:", format(stats_report$summary$unique_genomic_positions, big.mark = ",")),
        paste("Number of genes involved:", format(stats_report$summary$unique_genes, big.mark = ",")),
        paste("Number of cancer types:", stats_report$summary$unique_cancer_types),
        paste("Number of studies:", stats_report$summary$unique_studies),
        "",
        "=== Chromosome Distribution (Top 15) ===",
        capture.output(print(head(stats_report$chromosome_distribution, 15), row.names = FALSE)),
        "",
        "=== Cancer Type Distribution (Top 15) ===",
        capture.output(print(head(stats_report$cancer_type_distribution, 15), row.names = FALSE)),
        "",
        "=== Highly Mutated Genes (Top 20) ===",
        capture.output(print(head(stats_report$gene_distribution, 20), row.names = FALSE))
    )
    
    writeLines(summary_text, "final_summary_report.txt")
    print("Final summary report saved to: final_summary_report.txt")
    
    # Display key statistics
    print("")
    print("=== Key Statistics ===")
    print(paste("Most common chromosome:", stats_report$chromosome_distribution$chr[1], 
                "(mutations:", stats_report$chromosome_distribution$n[1], ")"))
    print(paste("Most common cancer type:", stats_report$cancer_type_distribution$cancer_type[1], 
                "(mutations:", stats_report$cancer_type_distribution$n[1], ")"))
    print(paste("Most common mutated gene:", stats_report$gene_distribution$hugoGeneSymbol[1], 
                "(mutations:", stats_report$gene_distribution$n[1], ")"))
}

# 11. Data validation function ------------------------------------------------------------
validate_results <- function(result) {
    if (is.null(result)) {
        print("No results to validate")
        return(FALSE)
    }
    
    print("Validating results...")
    
    mutations <- result$mutations
    genomic_positions <- result$genomic_positions
    
    # Check basic integrity
    checks <- list(
        has_mutations = !is.null(mutations) && nrow(mutations) > 0,
        has_genomic_positions = !is.null(genomic_positions) && nrow(genomic_positions) > 0,
        has_required_columns = all(c("chr", "startPosition", "endPosition", "hugoGeneSymbol") %in% names(genomic_positions)),
        no_na_positions = all(!is.na(genomic_positions$chr) & !is.na(genomic_positions$startPosition))
    )
    
    print("Validation results:")
    for (check_name in names(checks)) {
        status <- ifelse(checks[[check_name]], "✓", "✗")
        print(paste(" ", status, check_name))
    }
    
    # Display data preview
    if (checks$has_genomic_positions) {
        print("Genomic position data preview:")
        print(head(genomic_positions[, c("chr", "startPosition", "endPosition", "hugoGeneSymbol", "cancer_type")]))
    }
    
    return(all(unlist(checks)))
}

# 12. Function to list generated files ------------------------------------------------------
list_generated_files <- function() {
    print("List of generated files:")
    
    data_files <- c(
        "all_cancer_mutations_complete.csv",
        "all_cancer_mutations_genomic_positions.csv",
        "all_cancer_mutations_genome_browser.bed",
        "all_cancer_mutations_cancer_summary.csv",
        "all_cancer_mutations_gene_summary.csv",
        "mutation_statistics_summary.csv",
        "chromosome_distribution.csv",
        "cancer_type_distribution.csv",
        "gene_distribution.csv",
        "final_summary_report.txt"
    )
    
    # Check if mutation type distribution file exists
    if (file.exists("mutation_type_distribution.csv")) {
        data_files <- c(data_files, "mutation_type_distribution.csv")
    }
    
    for (file in data_files) {
        if (file.exists(file)) {
            file_info <- file.info(file)
            file_size <- round(file_info$size/1024/1024, 2)
            print(paste("✓", file, "-", file_size, "MB"))
        }
    }
    
    # Check cancer type files directory
    if (dir.exists("cancer_type_files")) {
        cancer_files <- list.files("cancer_type_files", pattern = "\\.csv$")
        print(paste("Cancer type files directory contains", length(cancer_files), "files"))
    }
}

# 13. Execute main program --------------------------------------------------------------

# Option 1: Test run (process first 10 studies)
#print("Starting test run (first 10 studies)...")
#test_result <- main(max_studies = 10, batch_size = 3)

# Validate test results
#if (exists("test_result")) {
#    test_validation <- validate_results(test_result)
 #   if (test_validation) {
 #       print("Test run validation passed!")
 #       list_generated_files()
#        
        # If test successful, ask whether to process all studies
 #       cat("\nTest run successful! Continue processing all remaining studies? (y/n): ")
  #      user_input <- readline()
 #       
#        if (tolower(user_input) == "y") {
#            print("Starting to process all studies...")
#            final_result <- main(max_studies = NULL, batch_size = 5)
#            
#            if (exists("final_result")) {
#                final_validation <- validate_results(final_result)
#                if (final_validation) {
#                    print("All studies processed and validation passed!")
#                    list_generated_files()
#                } else {
#                    print("All studies processed but validation failed!")
#                }
#            }
#        } else {
#            print("Only test run completed, not processing all studies")
#        }
#    } else {
#        print("Test run validation failed!")
#    }
#}

# Option 2: Directly process all studies (uncomment the following lines)
print("Starting to process all studies...")
final_result <- main(max_studies = NULL, batch_size = 5)

if (exists("final_result")) {
    final_validation <- validate_results(final_result)
    if (final_validation) {
        print("All studies processed and validation passed!")
        list_generated_files()
    } else {
        print("All studies processed but validation failed!")
    }
}

print("=== Program execution completed ===")

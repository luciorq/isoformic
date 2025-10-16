create_salmon_env <- function(env_name = "salmon-env", force = FALSE) {
  rlang::check_installed("condathis")
  condathis::create_env(
    package = "bioconda::salmon>=1.10.3",
    env_name = "salmon-env",
    additional_channels = "bioconda",
    verbose = "output",
    overwrite = force
  )
}

check_env_installed <- function(env_name = "salmon-env") {
  rlang::check_installed("condathis")
  rlang::check_required(env_name)
  if (isFALSE(condathis::env_exists(env_name))) {
    rlang::try_fetch(
      expr = {
        create_salmon_env(force = FALSE, env_name = env_name)
      },
      error = function(e) {
        cli::cli_abort(
          c(
            "The conda environment '{env_name}' could not be created or installed.",
            "i" = "Please try running {.code create_salmon_env()} manually.",
            "i" = "Original error message: {e$message}"
          ),
          class = "isoformic_salmon_env_create_error"
        )
      }
    )
  }
  return(invisible(TRUE))
}

salmon_index <- function(
  fasta_path,
  index_path = "salmon_index",
  kmer_len = 31,
  num_threads = 8,
  env_name = "salmon-env"
) {
  rlang::check_installed("condathis")
  rlang::check_required(fasta_path)
  check_env_installed(env_name)
  condathis::run(
    "salmon",
    "index",
    "--transcripts",
    fasta_path,
    "--index",
    index_path,
    "--type",
    "quasi",
    "--kmerLen",
    kmer_len,
    "--threads",
    num_threads,
    env_name = "salmon-env",
    verbose = "output"
  )
}

salmon_quant <- function(
  input_r1,
  input_r2 = NULL,
  index_path = "salmon_index",
  output_dir = "quant_output",
  num_threads = 8,
  num_gibbs = 100,
  env_name = "salmon-env"
) {
  rlang::check_installed("condathis")
  check_env_installed(env_name)
  rlang::check_required(input_r1)
  if (isFALSE(fs::dir_exists(index_path))) {
    cli::cli_abort(
      c(
        "The Salmon index directory '{index_path}' does not exist.",
        "i" = "Please run salmon_index() first to create the index."
      ),
      class = "isoformic_salmon_missing_index_error"
    )
  }

  condathis::run(
    "salmon",
    "quant",
    "--libType",
    "A",
    "--index",
    index_path,
    "--mates1",
    input_r1,
    "--mates2",
    input_r2,
    "--output",
    output_dir,
    "--threads",
    num_threads,
    "--validateMappings",
    "--d",
    "--posBias",
    "--seqBias",
    "--gcBias",
    "--numGibbsSamples",
    num_gibbs,
    env_name = "salmon-env",
    verbose = "output"
  )
}

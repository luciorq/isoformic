withr::local_options(
  .new = list(
    warnPartialMatchArgs = TRUE,
    warnPartialMatchAttr = TRUE,
    warnPartialMatchDollar = TRUE
  ),
  .local_envir = testthat::teardown_env()
)

tmp_home_path <- withr::local_tempdir(
  pattern = "tmp-home",
  .local_envir = testthat::teardown_env()
)
tmp_data_path <- withr::local_tempdir(
  pattern = "tmp-data",
  .local_envir = testthat::teardown_env()
)
tmp_cache_path <- withr::local_tempdir(
  pattern = "tmp-cache",
  .local_envir = testthat::teardown_env()
)

tmp_wd_path <- fs::path(tmp_home_path, "wd")

if (isFALSE(fs::dir_exists(tmp_home_path))) {
  fs::dir_create(tmp_home_path)
}
if (isFALSE(fs::dir_exists(tmp_data_path))) {
  fs::dir_create(tmp_data_path)
}
if (isFALSE(fs::dir_exists(tmp_cache_path))) {
  fs::dir_create(tmp_cache_path)
}
if (isFALSE(fs::dir_exists(tmp_wd_path))) {
  fs::dir_create(tmp_wd_path)
}

withr::local_dir(
  new = tmp_wd_path,
  .local_envir = testthat::teardown_env()
)

withr::local_envvar(
  .new = list(
    `HOME` = tmp_home_path,
    `USERPROFILE` = tmp_home_path,
    `LOCALAPPDATA` = tmp_data_path,
    `APPDATA` = tmp_data_path,
    `R_USER_DATA_DIR` = tmp_data_path,
    `XDG_DATA_HOME` = tmp_data_path,
    `R_USER_CACHE_DIR` = tmp_cache_path,
    `XDG_CACHE_HOME` = tmp_cache_path
  ),
  .local_envir = testthat::teardown_env()
)

base::rm(tmp_home_path)
base::rm(tmp_data_path)
base::rm(tmp_cache_path)
base::rm(tmp_wd_path)

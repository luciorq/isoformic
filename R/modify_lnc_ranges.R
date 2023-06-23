#' Edit the lncrna start and end table
#' Table of lncRNA starts and ends plus a value to add
#' @export
modify_lnc_ranges <- function(lncRNA_start_end, value_to_add) {
lncRNAisos_start_and_end_1mb <-lncRNA_start_end %>%
  mutate(added_value =value_to_add)

lncRNAisos_start_and_end_1mb_s <- lncRNAisos_start_and_end_1mb %>% 
  dplyr::select(-5) %>%
  mutate(new_start = c(TXSTART-added_value)) %>% 
  mutate(TXEND = lncRNAisos_start_and_end_1mb$TXEND)

lncRNAisos_start_and_end_1mb_e <- lncRNAisos_start_and_end_1mb_s %>% 
  dplyr::select(-6) %>%
  mutate(new_end = c(TXEND+added_value)) %>%
  mutate(new_start = lncRNAisos_start_and_end_1mb_s$new_start)

selected_lnc_newdist <- lncRNAisos_start_and_end_1mb_e %>% 
  dplyr::select(1:3, 7,8)
selected_lnc_newdist
}
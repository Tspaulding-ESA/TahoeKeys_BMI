########## FUNCTIONS #######################################################

`%notin%` <- Negate(`%in%`)

fix_names <- function(data){
  new_names = gsub("-|\\s+", "_", tolower(names(data)))
  setNames(data, new_names)
}

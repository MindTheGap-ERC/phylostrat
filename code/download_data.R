files = osfr::osf_retrieve_node("https://osf.io/zbpwa/") |>
  osfr::osf_ls_files() |>
  osfr::osf_download(path = "data/osf", conflicts = "overwrite")

# render all files
quarto::quarto_render("docs/graphs-tables.qmd", output_format = "all")
quarto::quarto_render("docs/supplemental-materials.qmd", output_format = "all")

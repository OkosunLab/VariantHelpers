filter_vars <- function(VCF,
                        biotype = c("protein_coding"),
                        impact = c("HIGH", "MODERATE"),
                        existing = TRUE,
                        PopulationFreq = 0.01,
                        AF = 0,
                        DP = 0,
                        ADP = 0
) {
    VCF <- mutate(VCF,
           Filter.biotype = BIOTYPE %in% biotype,
           Filter.impact = IMPACT %in% impact,
           Filter.existing = (! grepl("rs", Existing_variation) | grepl("COS", ID) ),
           Filter.population = (
               (AF_ALL <= PopulationFreq | is.na(AF_ALL)) &
                   (gnomADe_AF <= PopulationFreq | is.na(gnomADe_AF)) &
                   (gnomADg_AF <= PopulationFreq | is.na(gnomADg_AF))
               ),
           Filter.vaf = AF >= AF,
           Filter.depth = DP >= DP,
           Filter.alt.depth = ADP >= ADP,
           )
}

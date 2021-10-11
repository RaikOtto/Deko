aka3 = list(
    #Functionality   = c(
    #    ACTH = "white",
    #    PPoma = "BLACK",
    #    Unknown = "gray",
    #    Gastrinoma = "brown",
    #    Glucagonoma = "blue",
    #    Insulinoma = "darkgreen",
    #    "Non-functional" = "white",
    #    Somatostatinoma = "purple",
    #    VIPoma = "cyan"
    #),
    Evo_up_vec = c("High" = "red", "Medium" = "yellow","Low" = "green"),
    Functionality   = c(
        Functional = "green",
        Gastrinoma = "brown",
        Glucagonoma = "blue",
        Insulinoma = "darkgreen",
        Unknown = "gray",
        "Non-Functional" = "white",
        Ppoma = "purple",
        VIPoma = "cyan",
        Unknown = "gray"
    ),
    NEC_NET = c(NEC= "red", NET = "blue", Ambiguous = "purple", Unknown = "gray"),
    NET_NEC_PCA = c(NEC= "darkred", NET = "blue"),
    NEC_NET_PCA = c(NEC= "darkred", NET = "blue"),
    Study = c(Charite = "#35A047", Scarpa = "#F28500","Master" = "#008080"),
    Histology_Primary   = c(
        Pancreatic = "black",
        Colorectal = "brown",
        Intestinal = "white",
        Small_intestinal = "Yellow",
        Hepatic = "orange",
        Gastric = "purple",
        Gastric_not_specified = "gray"),
    Ratio = c("low"="white","high"="darkred"),
    Alpha = c(low = "white", high = "blue"),
    Beta = c(low = "white", high = "darkgreen"),
    Gamma = c(low = "white", high = "orange"),
    Delta = c(low = "white", high = "Purple"),
    Ductal = c(low = "white", high = "Brown"),
    Acinar = c(low = "white", high = "Cyan"),
<<<<<<< HEAD
    Grading = c( G1 = "#CBDB34",G2 = "#FFC000", G3 = "#EE5124"),
    Senescence_UP =c( "high" = "black", "low" = "white"),
    "TIS down" = c("high" = "darkred", "low" = "white"),
    P_value = c("not_sig" = "black","sig" = "white"),
    pannet_cluster = c(Left = "black", Right = "white"),
    Evo_up = c("High" = "red", "Medium" = "yellow", "Low" = "green"),
    CDH5 = c("High" = "red", "Medium" = "yellow", "Low" = "green"),
    Cohort = c("Left" = "black", Right = "white"),
    "PAKT_response" = c("Non-Responder_pAKT" = "red","Responder_pAKT" = "green"),
    "HER2" = c("HER2-" = "black", "HER2+" = "white", Unknown = "gray" )
=======
    Grading = c( G1 = "white",G2 = "#FFC000", G3 = "#EE5124", "Unknown" = "gray"),
    Primary_Metastasis = c("Primary" = "white", Metastasis = "gray", "Local_recurrence" = "pink"),
    MEN1 = c("0" = "white", "1"  = "#4073FF")
>>>>>>> 65dc5a559049fd0108f40cac2499cad5d57503a7
)

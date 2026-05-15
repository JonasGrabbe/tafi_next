####################################################################################
meta_names <- c(
  "site_index", "sample_index", "site_name", "sum", "time_range", "geological_context",
  "number_of_individuals", "age", "gender", "accumulation_cause", "hominin_species"
)
####################################################################################
sections_order <- c("cranium","teeth","hyoid", "vertebrae", "thorax", "shoulder",
"arms", "hands", "pelvic", "legs", "foot")

####################################################################################
high_resolution <- c(
  "Zanclean", "Piacenzian",
  "Gelasian", "Calabrian", "Chibanian", "Tarantian",
  "Oldowan", "Acheulean", "Mousterian", "Aterian",
  "Aurignacian", "Gravettian", "Solutrean", "Magdalenian",
  "Natufian", "Kebaran", "Azilian", "Maglemosian",
  "Pre-Pottery Neolithic A", "Pre-Pottery Neolithic B", "Pottery Neolithic",
  "Chalcolithic",  # Added this line
  "Early Bronze Age", "Middle Bronze Age", "Late Bronze Age",
  "Early Iron Age", "Middle Iron Age", "Late Iron Age",
  "Medieval Period", "Early Modern Period", "Late Modern Period"
)

low_resolution <- c(
  rep("Pliocene", 2),
  rep("Lower Pleistocene", 2), "Middle Pleistocene", "Upper Pleistocene",
  rep("Lower Paleolithic", 2), rep("Middle Paleolithic", 2),
  rep("Upper Paleolithic", 4),
  rep("Epipaleolithic/Mesolithic", 4),
  rep("Neolithic", 3),
  "Chalcolithic",
  rep("Bronze Age", 3),
  rep("Iron Age", 3),
  rep("Modern", 3)
)

# Define time ranges for each period (in years before present)
time_ranges <- list(

  #"Lower Paleolithic" = c(3400000, 300000),
  "Middle Paleolithic" = c(300000, 50000),
  "Upper Paleolithic" = c(50000, 11700),
  "Epipaleolithic/Mesolithic" = c(15000, 5000),
  "Neolithic" = c(12000, 5000),
  "Chalcolithic" = c(7000, 3300),
  "Bronze Age" = c(5300, 1200),
  "Iron Age" = c(3200, 500),
  "Modern" = c(1500, 0),


  "Pliocene" = c(5333000, 2580000),
  "Lower Pleistocene" = c(2580000, 781000),
  "Middle Pleistocene" = c(781000, 126000),
  #"Middle Pleistocene (before 0.3 MA)" = c(781000, 126000),
  "Upper Pleistocene" = c(126000, 11700)
)
# Define high-resolution time ranges for each period (in years before present)
high_resolution_time_ranges <- list(
  "Zanclean" = c(5333000, 3600000),      # 5.333-3.6 MA
  "Piacenzian" = c(3600000, 2580000),    # 3.6-2.58 MA
  "Gelasian" = c(2580000, 1800000),      # 2.58-1.8 MA
  "Calabrian" = c(1800000, 773000),      # 1.8-0.773 MA
  "Chibanian" = c(773000, 126000),       # 773-126 ka
  "Tarantian" = c(126000, 11700),        # 126-11.7 ka
  "Oldowan" = c(2600000, 1700000),       # 2.6-1.7 MA
  "Acheulean" = c(1700000, 300000),      # 1.7-0.3 MA
  "Mousterian" = c(300000, 30000),       # 300-30 ka
  "Aterian" = c(145000, 30000),          # 145-30 ka
  "Aurignacian" = c(43000, 26000),       # 43-26 ka
  "Gravettian" = c(26000, 21000),        # 26-21 ka
  "Solutrean" = c(21000, 17000),         # 21-17 ka
  "Magdalenian" = c(17000, 11700),       # 17-11.7 ka
  "Natufian" = c(15000, 12500),          # 15-12.5 ka
  "Kebaran" = c(18000, 12500),           # 18-12.5 ka
  "Azilian" = c(11000, 9000),            # 11-9 ka
  "Maglemosian" = c(9000, 5800),         # 9-5.8 ka
  "Pre-Pottery Neolithic A" = c(10000, 9000),  # 10-9 ka
  "Pre-Pottery Neolithic B" = c(9000, 8500),   # 9-8.5 ka
  "Pottery Neolithic" = c(8500, 7000),         # 8.5-7 ka
  "Early Bronze Age" = c(5300 + 1950, 2100 + 1950),  # 3300-2100 BC
  "Middle Bronze Age" = c(2100 + 1950, 1550 + 1950), # 2100-1550 BC
  "Late Bronze Age" = c(1550 + 1950, 1200 + 1950),   # 1550-1200 BC
  "Early Iron Age" = c(1200 + 1950, 800 + 1950),     # 1200-800 BC
  "Middle Iron Age" = c(800 + 1950, 500 + 1950),     # 800-500 BC
  "Late Iron Age" = c(500 + 1950, 0),               # 500-0 BC
  "Medieval Period" = c(500, 1500),                 # 500-1500 AD
  "Early Modern Period" = c(1500, 1800),            # 1500-1800 AD
  "Late Modern Period" = c(1800, 0)                 # 1800 AD-Present
)
####################################################################################
hominin_species <- c(
  "Australopithecus anamensis",   # 4.2 - 3.9 million years ago
  "Australopithecus afarensis",   # 3.9 - 2.9 million years ago
  "Australopithecus africanus",   # 3.0 - 2.1 million years ago
  "Paranthropus aethiopicus",     # 2.7 - 2.3 million years ago
  "Homo habilis",                 # 2.8 - 1.5 million years ago
  "Homo rudolfensis",             # 2.4 - 1.9 million years ago
  "Paranthropus boisei",          # 2.3 - 1.2 million years ago
  "Paranthropus robustus",        # 2.0 - 1.2 million years ago
  "Homo ergaster",                # 1.9 - 1.4 million years ago
  "Homo erectus",                 # 1.9 million - 110,000 years ago
  "Homo antecessor",              # 1.2 million - 800,000 years ago
  "Homo heidelbergensis",         # 600,000 - 200,000 years ago
  "Neanderthals (Homo neanderthalensis)",  # 400,000 - 40,000 years ago
  "Denisovans",                   # 500,000 - 30,000 years ago
  "Homo naledi",                  # 335,000 - 236,000 years ago
  "Homo sapiens"                  # 300,000 years ago - present
)
####################################################################################
geological_context <- c(
  "Cave",
  "Sepulchral Cave",  # Included as a specific type of cave
  "Rockshelter",
  "Open air",
  "Surface disposal"
)
####################################################################################
age_category <- c(
  "Infant (0-2 years)",
  "Child (3-7 years)",
  "Adolescent (8-16 years)",
  "Subadult",  # General category if specific age is not available
  "Adult (17+ years)"
)
##########################################

##########################################
gender <- c(
  "Male",
  "Female",
  "Unknown" # Mixed groups or indeterminate
)
####################################################################################
accumulation_cause <- c(
  "Cannibalism",
  "Burial",
  "Collective burial",
  "Mass graves",  # Specific large-scale collective burials
  "Secondary position",  # Remains moved post-mortem
  "Surface disposal",  # Remains left exposed
  "Disposal",  # General category for any non-ceremonial discarding
  "Scavenged",  # Altered by animals post-mortem
  "Intentional accumulation",  # Human-driven collection
  "Funerary assemblage",  # Group of burial goods and remains
  "Cemetery",  # Organized burial area
  "Cuevas Sepulcrales",  # Sepulchral caves, specific burial sites
  "Transport or accumulation from carnivores",  # Non-human driven accumulation
  "Possible burial"  # Uncertain if burial context
)
####################################################################################
# Define the hierarchy for each section with bone counts
skeleton_hierarchy <- list(

  cranium = list(
    "Level 1" = list(
      cranium = 1
    ),

    "Level 2" = list(
      neurocranium_bone = 1,
      splanchnocranium  = 1
    ),

    "Level 3" = list(
      neurocranium_bone_indet = 1,
      frontal_bone  = 1,
      parietal_bone = 2,
      temporal      = 2,
      occipital     = 1,
      sphenoid      = 1,
      ethmoid       = 1,

      splanchnocranium_indet = 1,
      maxilla               = 2,
      palatine_bone         = 2,
      zygomatic_bone        = 2,
      nasal_bone            = 2,
      lacrimal_bone         = 2,
      vomer_bone            = 1,
      inferior_nasal_concha = 2,

      ossicles_of_the_middle_ear = 6,
      mandible = 1
    ),

    "Level 4" = list(
      # carried down (min=3 max=5 in levels_down (1).xlsx)
      neurocranium_bone_indet = 1,
      

      frontal_bone = 1,

      parietal_bone_indet = 2,
      parietal_bone_right = 1,
      parietal_bone_left  = 1,

      temporal_bone_indet = 2,
      temporal_right      = 1,
      temporal_left       = 1,

      occipital = 1,
      sphenoid  = 1,
      ethmoid   = 1,

      splanchnocranium_indet  = 1,
      maxilla       = 2,
      palatine_bone = 2,

      zygomatic_bone_indet = 2,
      zygomatic_bone_right = 1,
      zygomatic_bone_left  = 1,

      nasal_bone    = 2,
      lacrimal_bone = 2,
      vomer_bone    = 1,
      inferior_nasal_concha = 2,

      ossicles_of_the_middle_ear_indet = 6,
      malleus = 2,
      incus   = 2,
      stapes  = 2,

      mandible = 1
    ),

    "Level 5" = list(
      # carried down (min=3 max=5 in levels_down (1).xlsx)
      neurocranium_bone_indet = 1,
      

      frontal_bone = 1,

      parietal_bone_indet = 2,
      parietal_bone_right = 1,
      parietal_bone_left  = 1,

      temporal_bone_indet = 2,
      temporal_right      = 1,
      temporal_left       = 1,

      occipital = 1,
      sphenoid  = 1,
      ethmoid   = 1,

      splanchnocranium_indet  = 1,
      maxilla       = 2,
      palatine_bone = 2,

      zygomatic_bone_indet = 2,
      zygomatic_bone_right = 1,
      zygomatic_bone_left  = 1,

      nasal_bone    = 2,
      lacrimal_bone = 2,
      vomer_bone    = 1,
      inferior_nasal_concha = 2,

      ossicles_of_the_middle_ear_indet = 6,

      malleus_indet = 2,
      malleus_right = 1,
      malleus_left  = 1,

      incus_indet = 2,
      incus_right = 1,
      incus_left  = 1,

      stapes_indet = 2,
      stapes_right = 1,
      stapes_left  = 1,

      mandible = 1
    )
  ),

  teeth = list(
    "Level 1" = list(
      teeth = 1
    ),

    "Level 2" = list(
      #deciduous_indet = 1, 
      deciduous = 20,
      permanent = 32
    ),

    "Level 3" = list(
      #deciduous_indet = 0,
      #permanent_indet = 0,

      di1 = 4,
      di2 = 4,
      dc  = 4,
      dm1 = 4,
      dm2 = 4,

      i1  = 4,
      i2  = 4,
      c   = 4,
      pm3 = 4,
      pm4 = 4,
      m1  = 4,
      m2  = 4,
      m3  = 4
    ),

    "Level 4" = list(
      #deciduous_indet = 0,
      #permanent_indet = 0,

      dl1_indet     = 4,
      di1_superior  = 2,
      dl1_inferior  = 2,

      dl2_indet     = 4,
      di2_superior  = 2,
      dl2_inferior  = 2,

      dc_indet      = 4,
      dc_superior   = 2,
      dc_inferior   = 2,

      dm1_indet     = 4,
      dm1_superior  = 2,
      dm1_inferior  = 2,

      dm2_indet     = 4,
      dm2_superior  = 2,
      dm2_inferior  = 2,

      i1_indet      = 4,
      i1_superior   = 2,
      i1_inferior   = 2,

      i2_indet      = 4,
      i2_superior   = 2,
      i2_inferior   = 2,

      c_indet       = 4,
      c_superior    = 2,
      c_inferior    = 2,

      pm3_indet     = 4,
      pm3_superior  = 2,
      pm3_inferior  = 2,

      pm4_indet     = 4,
      pm4_superior  = 2,
      pm4_inferior  = 2,

      m1_indet      = 4,
      m1_superior   = 2,
      m1_inferior   = 2,

      m2_indet      = 4,
      m2_superior   = 2,
      m2_inferior   = 2,

      m3_indet      = 4,
      m3_superior   = 2,
      m3_inferior   = 2
    ),

    "Level 5" = list(
      #deciduous_indet = 0,
      #permanent_indet = 0,

      dl1_indet = 4,
      di1_superior_indet = 2,
      di1_superior_right = 1,
      di1_superior_left  = 1,
      dl1_inferior_indet = 2,
      dl1_inferior_right = 1,
      dl1_inferior_left  = 1,

      dl2_indet = 4,
      di2_superior_indet = 2,
      di2_superior_right = 1,
      di2_superior_left  = 1,
      di2_inferior_indet = 2,
      di2_inferior_right = 1,
      di2_inferior_left  = 1,

      dc_indet = 4,
      dc_superior_indet = 2,
      dc_superior_right = 1,
      dc_superior_left  = 1,
      dc_inferior_indet = 2,
      dc_inferior_right = 1,
      dc_inferior_left  = 1,

      dm1_indet = 4,
      dm1_superior_indet = 2,
      dm1_superior_right = 1,
      dm1_superior_left  = 1,
      dm1_inferior_indet = 2,
      dm1_inferior_right = 1,
      dm1_inferior_left  = 1,

      dm2_indet = 4,
      dm2_superior_indet = 2,
      dm2_superior_right = 1,
      dm2_superior_left  = 1,
      dm2_inferior_indet = 2,
      dm2_inferior_right = 1,
      dm2_inferior_left  = 1,

      i1_indet = 4,
      i1_superior_indet = 2,
      i1_superior_right = 1,
      i1_superior_left  = 1,
      i1_inferior_indet = 2,
      i1_inferior_right = 1,
      i1_inferior_left  = 1,

      i2_indet = 4,
      i2_superior_indet = 2,
      i2_superior_right = 1,
      i2_superior_left  = 1,
      i2_inferior_indet = 2,
      i2_inferior_right = 1,
      i2_inferior_left  = 1,

      c_indet = 4,
      c_superior_indet = 2,
      c_superior_right = 1,
      c_superior_left  = 1,
      c_inferior_indet = 2,
      c_inferior_right = 1,
      c_inferior_left  = 1,

      pm3_indet = 4,
      pm3_superior_indet = 2,
      pm3_superior_right = 1,
      pm3_superior_left  = 1,
      pm3_inferior_indet = 2,
      pm3_inferior_right = 1,
      pm3_inferior_left  = 1,

      pm4_indet = 4,
      pm4_superior_indet = 2,
      pm4_superior_right = 1,
      pm4_superior_left  = 1,
      pm4_inferior_indet = 2,
      pm4_inferior_right = 1,
      pm4_inferior_left  = 1,

      m1_indet = 4,
      m1_superior_indet = 2,
      m1_superior_right = 1,
      m1_superior_left  = 1,
      m1_inferior_indet = 2,
      m1_inferior_right = 1,
      m1_inferior_left  = 1,

      m2_indet = 4,
      m2_superior_indet = 2,
      m2_superior_right = 1,
      m2_superior_left  = 1,
      m2_inferior_indet = 2,
      m2_inferior_right = 1,
      m2_inferior_left  = 1,

      m3_indet = 4,
      m3_superior_indet = 2,
      m3_superior_right = 1,
      m3_superior_left  = 1,
      m3_inferior_indet = 2,
      m3_inferior_right = 1,
      m3_inferior_left  = 1
    )
  ),

  hyoid = list(
    "Level 1" = list(
      hyoid = 1
    )
  ),

  vertebrae = list(
    "Level 1" = list(
      vertebrae = 24
    ),

    "Level 2" = list(
      vertebrae_indet = 24,
      cervical = 7,
      thoracic = 12,
      lumbar   = 5
    ),

    "Level 3" = list(
      vertebrae_indet = 24,

      cervical_indet = 7,
      c1_atlas = 1,
      c2_axis  = 1,
      c3 = 1,
      c4 = 1,
      c5 = 1,
      c6 = 1,
      c7 = 1,

      thoracic_indet = 12,
      t1  = 1,
      t2  = 1,
      t3  = 1,
      t4  = 1,
      t5  = 1,
      t6  = 1,
      t7  = 1,
      t8  = 1,
      t9  = 1,
      t10 = 1,
      t11 = 1,
      t12 = 1,

      lumbar_indet = 5,
      l1 = 1,
      l2 = 1,
      l3 = 1,
      l4 = 1,
      l5 = 1
    )
  ),

  thorax = list(
    "Level 1" = list(
      thorax = 25
    ),

    "Level 2" = list(
      thorax_indet = 25,
      sternum = 1,
      ribs = 24
    ),

    "Level 3" = list(
      # carried down (min=2 max=4)
      thorax_indet = 25,
      sternum = 1,

      rib_indet     = 24,
      rib_fragments = 72,
      

      ribs_right = 12,
      ribs_left = 12,

      rib_1  = 2,
      rib_2  = 2,
      rib_3  = 2,
      rib_4  = 2,
      rib_5  = 2,
      rib_6  = 2,
      rib_7  = 2,
      rib_8 = 2,
      rib_9  = 2,
      rib_10 = 2,
      rib_11 = 2,
      rib_12 = 2
    ),

    "Level 4" = list(
      # carried down (min=2 max=4)
      thorax_indet = 25,
      sternum = 1,

      # carried down (min=3 max=4)
      rib_indet     = 24,
      rib_fragments = 72,
      
      ribs_right = 12,
      ribs_left = 12,

      #rib_1_indet  = 0, 
      rib_1_right  = 1,
      rib_1_left   = 1,
      
      #rib_2_indet = 0,
      rib_2_right  = 1,
      rib_2_left   = 1,

      #rib_3_indet = 0,
      rib_3_right  = 1,
      rib_3_left   = 1,

      #rib_4_indet = 0,
      rib_4_right  = 1,
      rib_4_left   = 1,

      #rib_5_indet = 0,
      rib_5_right  = 1,
      rib_5_left   = 1,

      #rib_6_indet = 0,
      rib_6_right  = 1,
      rib_6_left   = 1,

      #rib_7_indet = 0,
      rib_7_right  = 1,
      rib_7_left   = 1,

      #rib_8_indet = 0,
      rib_8_right  = 1,
      rib_8_left   = 1,

      #rib_9_indet = 0,
      rib_9_right  = 1,
      rib_9_left   = 1,

      #rib_10_indet = 0,
      rib_10_right = 1,
      rib_10_left  = 1,

      #rib_11_indet = 0,
      rib_11_right = 1,
      rib_11_left  = 1,

      #rib_12_indet = 0,
      rib_12_right = 1,
      rib_12_left  = 1
    )
  ),

  shoulder = list(
    "Level 1" = list(
      shoulder = 4
    ),

    "Level 2" = list(
      shoulder_indet = 4,
      clavicle = 2,
      scapula  = 2
    ),

    "Level 3" = list(
      # carried down (min=2 max=3)
      shoulder_indet = 4,

      clavicle_indet = 2,
      clavicle_right = 1,
      clavicle_left  = 1,

      scapula_indet = 2,
      scapula_right = 1,
      scapula_left  = 1
    )
  ),

  arms = list(
    "Level 1" = list(
      arms = 6
    ),

    "Level 2" = list(
      arms_indet = 6,
      humerus = 2,
      radius  = 2,
      ulna    = 2
    ),

    "Level 3" = list(
      # carried down (min=2 max=3)
      arms_indet = 6,

      humerus_indet = 2,
      humerus_right = 1,
      humerus_left  = 1,

      radius_indet = 2,
      radius_right = 1,
      radius_left  = 1,

      ulna_indet = 2,
      ulna_right = 1,
      ulna_left  = 1
    )
  ),

  hands = list(
    "Level 1" = list(
      hands = 54
    ),

    "Level 2" = list(
      hands_indet  = 54,
      carpal       = 16,
      metacarpal   = 10,
      hand_phalanx = 28
    ),

    "Level 3" = list(
      hands_indet = 54,

      carpal_indet = 16,
      carpal_right = 8,
      carpal_left  = 8,

      metacarpal_indet = 10,
      metacarpal_right = 5,
      metacarpal_left  = 5,

      hand_phalanx_indet = 28,
      hand_phalanx_proximal = 10,
      hand_phalanx_middle   = 8,
      hand_phalanx_distal      = 10
    ),

    "Level 4" = list(
      # carried down (min<=3 max=5)
      hands_indet = 54,
      carpal_indet = 16,
      metacarpal_indet = 10,
      hand_phalanx_indet = 28,

      carpal_right_indet = 8,
      scaphoid_right   = 1,
      lunate_right     = 1,
      triquetrum_right = 1,
      pisiform_right   = 1,
      trapezium_right  = 1,
      trapezoid_right  = 1,
      capitate_right   = 1,
      hamate_right     = 1,

      carpal_left_indet = 8,
      scaphoid_left   = 1,
      lunate_left     = 1,
      triquetrum_left = 1,
      pisiform_left   = 1,
      trapezium_left  = 1,
      trapezoid_left  = 1,
      capitate_left   = 1,
      hamate_left     = 1,

      metacarpal_right_indet = 5,
      metacarpal_1_right = 1,
      metacarpal_2_right = 1,
      metacarpal_3_right = 1,
      metacarpal_4_right = 1,
      metacarpal_5_right = 1,

      metacarpal_left_indet = 5,
      metacarpal_1_left = 1,
      metacarpal_2_left = 1,
      metacarpal_3_left = 1,
      metacarpal_4_left = 1,
      metacarpal_5_left = 1,

      hand_phalanx_proximal_indet = 10,
      hand_phalanx_right = 5,
      hand_phalanx_left  = 5,

      hand_phalanx_middle_indet = 8,
      hand_middle_right = 4,
      hand_middle_left  = 4,

      hand_phalanx_distal_indet = 10,
      hand_phalanx_distal_right = 5,
      hand_distal_left          = 5
    ),

    "Level 5" = list(
      # carry down Level 4 (min<=4 max=5) + add Level 5 specifics
      hands_indet = 54,
      carpal_indet = 16,
      metacarpal_indet = 10,
      hand_phalanx_indet = 28,

      carpal_right_indet = 8,
      scaphoid_right   = 1,
      lunate_right     = 1,
      triquetrum_right = 1,
      pisiform_right   = 1,
      trapezium_right  = 1,
      trapezoid_right  = 1,
      capitate_right   = 1,
      hamate_right     = 1,

      carpal_left_indet = 8,
      scaphoid_left   = 1,
      lunate_left     = 1,
      triquetrum_left = 1,
      pisiform_left   = 1,
      trapezium_left  = 1,
      trapezoid_left  = 1,
      capitate_left   = 1,
      hamate_left     = 1,

      metacarpal_right_indet = 5,
      metacarpal_1_right = 1,
      metacarpal_2_right = 1,
      metacarpal_3_right = 1,
      metacarpal_4_right = 1,
      metacarpal_5_right = 1,

      metacarpal_left_indet = 5,
      metacarpal_1_left = 1,
      metacarpal_2_left = 1,
      metacarpal_3_left = 1,
      metacarpal_4_left = 1,
      metacarpal_5_left = 1,

      hand_phalanx_proximal_indet = 10,
      hand_phalanx_proximal_right_indet = 5,
      hand_phalanx_proximal_left_indet  = 5,
      hand_phalanx_proximal_1_right = 1,
      hand_phalanx_proximal_2_right = 1,
      hand_phalanx_proximal_3_right = 1,
      hand_phalanx_proximal_4_right = 1,
      hand_phalanx_proximal_5_right = 1,
      hand_phalanx_proximal_1_left = 1,
      hand_phalanx_proximal_2_left = 1,
      hand_phalanx_proximal_3_left = 1,
      hand_phalanx_proximal_4_left = 1,
      hand_phalanx_proximal_5_left = 1,

      hand_phalanx_middle_indet = 8,
      hand_phalanx_middle_right_indet = 4,
      hand_phalanx_middle_left_indet  = 4,
      hand_phalanx_middle_2_right = 1,
      hand_phalanx_middle_3_right = 1,
      hand_phalanx_middle_4_right = 1,
      hand_phalanx_middle_5_right = 1,
      hand_phalanx_middle_2_left = 1,
      hand_phalanx_middle_3_left = 1,
      hand_phalanx_middle_4_left = 1,
      hand_phalanx_middle_5_left = 1,

      hand_phalanx_distal_indet = 10,
      hand_phalanx_distal_right_indet = 5,
      hand_phalanx_distal_left_indet  = 5,
      hand_phalanx_distal_right_1 = 1,
      hand_phalanx_distal_right_2 = 1,
      hand_phalanx_distal_right_3 = 1,
      hand_phalanx_distal_right_4 = 1,
      hand_phalanx_distal_right_5 = 1,
      hand_phalanx_distal_left_1 = 1,
      hand_phalanx_distal_left_2 = 1,
      hand_phalanx_distal_left_3 = 1,
      hand_phalanx_distal_left_4 = 1,
      hand_phalanx_distal_left_5 = 1
    )
  ),

  pelvic = list(
    "Level 1" = list(
      pelvic = 1
    ),

    "Level 2" = list(
      pelvic_indet = 1,
      sacrum  = 1,
      coccyx  = 1,
      os_coxae = 2
    ),

    "Level 3" = list(
      pelvic_indet = 1,

      sacrum_indet = 1,
      s1 = 1,
      s2 = 1,
      s3 = 1,
      s4 = 1,
      s5 = 1,

      coccyx_indet = 1,
      cx1 = 1,
      cx2 = 1,
      cx3 = 1,
      cx4 = 1,

      os_coxae_indet = 2,
      os_coxae_right = 1,
      os_coxae_left  = 1
    ),

    "Level 4" = list(
      pelvic_indet = 1,

      sacrum_indet = 1,
      s1 = 1,
      s2 = 1,
      s3 = 1,
      s4 = 1,
      s5 = 1,
      

      coccyx_indet = 1,
      cx1 = 1,
      cx2 = 1,
      cx3 = 1,
      cx4 = 1,
      

      os_coxae_indet = 2,
      os_coxae_right_indet = 1,
      ilium_right   = 1,
      ischium_right = 1,
      pubis_right   = 1,

      os_coxae_left_indet = 1,
      ilium_left   = 1,
      ischium_left = 1,
      pubis_left   = 1
    )
  ),
  legs = list(
    "Level 1" = list(
      legs = 8
    ),

    "Level 2" = list(
      legs_indet = 8,
      femur   = 2,
      tibia   = 2,
      fibula  = 2,
      patella = 2
    ),

    "Level 3" = list(
      legs_indet = 8,

      femur_indet = 2,
      femur_right = 1,
      femur_left  = 1,

      tibia_indet = 2,
      tibia_right = 1,
      tibia_left  = 1,

      fibula_indet = 2,
      fibula_right = 1,
      fibula_left  = 1,

      patella_indet = 2,
      patella_right = 1,
      patella_left  = 1
    )
  ),

  foot = list(
    "Level 1" = list(
      foot = 52
    ),

    "Level 2" = list(
      foot_indet    = 52,
      tarsal        = 14,
      metatarsal    = 10,
      foot_phalanx  = 28
    ),

    "Level 3" = list(
      foot_indet = 52,

      tarsal_indet = 14,
      tarsal_right = 7,
      tarsal_left  = 7,

      metatarsal_indet = 10,
      metatarsal_right = 5,
      metatarsal_left  = 5,

      foot_phalanx_indet = 28,
      foot_phalanx_proximal = 10,
      foot_phalanx_middle   = 8,
      foot_phalanx_distal      = 10
    ),

    "Level 4" = list(
      # carry-down (min<=4<=max)
      foot_indet = 52,

      tarsal_indet = 14,
      tarsal_right_indet = 7,
      talus_right        = 1,
      calcaneous_right   = 1,
      cuboid_right       = 1,
      navicular_right    = 1,
      lateral_cuneiform_right       = 1,
      intermediate_cuneiform_right  = 1,
      medial_cuneiform_right        = 1,

      tarsal_left_indet = 7,
      talus_left        = 1,
      calcaneous_left   = 1,
      cuboid_left       = 1,
      navicular_left    = 1,
      lateral_cuneiform_left       = 1,
      intermediate_cuneiform_left  = 1,
      medial_cuneiform_left        = 1,

      metatarsal_indet = 10,
      metatarsal_right_indet = 5,
      metatarsal_1_right = 1,
      metatarsal_2_right = 1,
      metatarsal_3_right = 1,
      metatarsal_4_right = 1,
      metatarsal_5_right = 1,

      metatarsal_left_indet = 5,
      metatarsal_1_left = 1,
      metatarsal_2_left = 1,
      metatarsal_3_left = 1,
      metatarsal_4_left = 1,
      metatarsal_5_left = 1,

      foot_phalanx_indet = 28,

      foot_phalanx_proximal_indet = 10,
      foot_proximal_right = 5,
      foot_phalanx_left   = 5,

      foot_phalanx_middle_indet = 8,
      foot_middle_right = 4,
      foot_middle_left  = 4,

      foot_phalanx_distal_indet = 10,
      foot_phalanx_distal_right = 5,
      foot_distal_left          = 5
    ),

    "Level 5" = list(
      foot_indet = 52,

      tarsal_indet = 14,
      tarsal_right_indet = 7,
      talus_right        = 1,
      calcaneous_right   = 1,
      cuboid_right       = 1,
      navicular_right    = 1,
      lateral_cuneiform_right       = 1,
      intermediate_cuneiform_right  = 1,
      medial_cuneiform_right        = 1,

      tarsal_left_indet = 7,
      talus_left        = 1,
      calcaneous_left   = 1,
      cuboid_left       = 1,
      navicular_left    = 1,
      lateral_cuneiform_left       = 1,
      intermediate_cuneiform_left  = 1,
      medial_cuneiform_left        = 1,

      metatarsal_indet = 10,
      metatarsal_right_indet = 5,
      metatarsal_1_right = 1,
      metatarsal_2_right = 1,
      metatarsal_3_right = 1,
      metatarsal_4_right = 1,
      metatarsal_5_right = 1,

      metatarsal_left_indet = 5,
      metatarsal_1_left = 1,
      metatarsal_2_left = 1,
      metatarsal_3_left = 1,
      metatarsal_4_left = 1,
      metatarsal_5_left = 1,

      foot_phalanx_indet = 28,

      foot_phalanx_proximal_indet = 10,
      foot_phalanx_proximal_right_indet = 5,
      foot_phalanx_proximal_left_indet  = 5,
     
      foot_phalanx_proximal_1_right = 1,
      foot_phalanx_proximal_2_right = 1,
      foot_phalanx_proximal_3_right = 1,
      foot_phalanx_proximal_4_right = 1,
      foot_phalanx_proximal_5_right = 1,
      foot_phalanx_proximal_1_left = 1,
      foot_phalanx_proximal_2_left = 1,
      foot_phalanx_proximal_3_left = 1,
      foot_phalanx_proximal_4_left = 1,
      foot_phalanx_proximal_5_left = 1,

      foot_phalanx_middle_indet = 8,
      foot_phalanx_middle_right_indet = 4,
      foot_phalanx_middle_left_indet  = 4,
      foot_phalanx_middle_2_right = 1,
      foot_phalanx_middle_3_right = 1,
      foot_phalanx_middle_4_right = 1,
      foot_phalanx_middle_5_right = 1,
      foot_phalanx_middle_2_left = 1,
      foot_phalanx_middle_3_left = 1,
      foot_phalanx_middle_4_left = 1,
      foot_phalanx_middle_5_left = 1,

      foot_phalanx_distal_indet = 10,
      foot_phalanx_distal_right_indet = 5,
      foot_phalanx_distal_left_indet  = 5,
      foot_phalanx_distal_right_1 = 1,
      foot_phalanx_distal_right_2 = 1,
      foot_phalanx_distal_right_3 = 1,
      foot_phalanx_distal_right_4 = 1,
      foot_phalanx_distal_right_5 = 1,
      foot_phalanx_distal_left_1 = 1,
      foot_phalanx_distal_left_2 = 1,
      foot_phalanx_distal_left_3 = 1,
      foot_phalanx_distal_left_4 = 1,
      foot_phalanx_distal_left_5 = 1
    )
  )

)
############################################################################################
ElementsOne <- c(
  
  "cranium", 
  "neurocranium_bone", "neurocranium_bone_indet", "frontal_bone", "parietal_bone", "parietal_bone_indet",
  "parietal_bone_right", "parietal_bone_left", "temporal", "temporal_bone_indet", "temporal_right", "temporal_left",
  "occipital", "sphenoid", "ethmoid", 
  "splanchnocranium", "splanchnocranium_indet", "maxilla", "palatine_bone",
  "zygomatic_bone", "zygomatic_bone_indet", "zygomatic_bone_right", "zygomatic_bone_left", "nasal_bone",
  "lacrimal_bone", "vomer_bone", "inferior_nasal_concha", "ossicles_of_the_middle_ear",
  "ossicles_of_the_middle_ear_indet", "malleus", "malleus_indet", "malleus_right", "malleus_left", "incus",
  "incus_indet", "incus_right", "incus_left", "stapes", "stapes_indet", "stapes_right", "stapes_left", 
  "mandible",
  
  
  "teeth", 
  "deciduous",
  #"deciduous_indet", 
  "di1", "dl1_indet", "di1_superior", "di1_superior_indet", "di1_superior_right",
  "di1_superior_left", "dl1_inferior", "dl1_inferior_indet", "dl1_inferior_right", "dl1_inferior_left", "di2",
  "dl2_indet", "di2_superior", "di2_superior_indet", "di2_superior_right", "di2_superior_left", "dl2_inferior",
  "di2_inferior_indet", "di2_inferior_right", "di2_inferior_left", "dc", "dc_indet", "dc_superior",
  "dc_superior_indet", "dc_superior_right", "dc_superior_left", "dc_inferior", "dc_inferior_indet",
  "dc_inferior_right", "dc_inferior_left", "dm1", "dm1_indet", "dm1_superior", "dm1_superior_indet",
  "dm1_superior_left", "dm1_superior_right", "dm1_inferior", "dm1_inferior_indet", "dm1_inferior_right",
  "dm1_inferior_left", "dm2", "dm2_indet", "dm2_superior", "dm2_superior_indet", "dm2_superior_right",
  "dm2_superior_left", "dm2_inferior", "dm2_inferior_indet", "dm2_inferior_right", "dm2_inferior_left", 
  
  "permanent", 
  #"permanent_indet",
  "i1", "i1_indet", "i1_superior", "i1_superior_indet", "i1_superior_right", "i1_superior_left", "i1_inferior",
  "i1_inferior_indet", "i1_inferior_right", "i1_inferior_left", "i2", "i2_indet", "i2_superior", "i2_superior_indet",
  "i2_superior_right", "i2_superior_left", "i2_inferior", "i2_inferior_indet", "i2_inferior_right", "i2_inferior_left",
  "c", "c_indet", "c_superior", "c_superior_indet", "c_superior_right", "c_superior_left", "c_inferior",
  "c_inferior_indet", "c_inferior_right", "c_inferior_left", "pm3", "pm3_indet", "pm3_superior", "pm3_superior_indet",
  "pm3_superior_right", "pm3_superior_left", "pm3_inferior", "pm3_inferior_indet", "pm3_inferior_right",
  "pm3_inferior_left", "pm4", "pm4_indet", "pm4_superior", "pm4_superior_indet", "pm4_superior_right",
  "pm4_superior_left", "pm4_inferior", "pm4_inferior_indet", "pm4_inferior_right", "pm4_inferior_left", "m1",
  "m1_indet", "m1_superior", "m1_superior_indet", "m1_superior_right", "m1_superior_left", "m1_inferior",
  "m1_inferior_indet", "m1_inferior_right", "m1_inferior_left", "m2", "m2_indet", "m2_superior", "m2_superior_indet",
  "m2_superior_right", "m2_superior_left", "m2_inferior", "m2_inferior_indet", "m2_inferior_right", "m2_inferior_left",
  "m3", "m3_indet", "m3_superior", "m3_superior_indet", "m3_superior_right", "m3_superior_left", "m3_inferior",
  "m3_inferior_indet", "m3_inferior_right", "m3_inferior_left", 
  
  "hyoid", 
  
  "vertebrae", "vertebrae_indet", "cervical",
  "cervical_indet", "c1_atlas", "c2_axis", "c3", "c4", "c5", "c6", "c7", "thoracic", "thoracic_indet", "t1", "t2",
  "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10", "t11", "t12", "lumbar", "lumbar_indet", "l1", "l2", "l3", "l4",
  "l5", 
  
  "thorax", 
  "thorax_indet", 
  "sternum", 
  "ribs", "rib_indet", "rib_fragments", "ribs_right", "ribs_left",
  "rib_1", "rib_1_right", "rib_1_left", "rib_2", "rib_2_right", "rib_2_left", "rib_3", "rib_3_right",
  "rib_3_left", "rib_4", "rib_4_right", "rib_4_left", "rib_5", "rib_5_right", "rib_5_left", "rib_6",
  "rib_6_right", "rib_6_left", "rib_7", "rib_7_right", "rib_7_left", "rib_8", "rib_8_right", "rib_8_left",
  "rib_9", "rib_9_right", "rib_9_left", "rib_10", "rib_10_right", "rib_10_left", "rib_11",
  "rib_11_right", "rib_11_left", "rib_12", "rib_12_right", "rib_12_left", 
  
  "shoulder", "shoulder_indet",
  "clavicle", "clavicle_indet", "clavicle_right", "clavicle_left", "scapula", "scapula_indet", "scapula_right",
  "scapula_left", 
  
  "arms", "arms_indet", 
  "humerus", "humerus_indet", "humerus_right", "humerus_left", "radius",
  "radius_indet", "radius_right", "radius_left", "ulna", "ulna_indet", "ulna_right", "ulna_left", 
  
  "hands",
  "hands_indet", 
  "carpal", "carpal_indet", "carpal_right", "carpal_right_indet", "trapezium_right", "trapezoid_right",
  "hamate_right", "triquetrum_right", "lunate_right", "pisiform_right", "capitate_right", "scaphoid_right",
  "carpal_left", "carpal_left_indet", "trapezium_left", "trapezoid_left", "hamate_left", "triquetrum_left",
  "lunate_left", "pisiform_left", "capitate_left", "scaphoid_left", 
  "metacarpal", "metacarpal_indet",
  "metacarpal_right", "metacarpal_right_indet", "metacarpal_1_right", "metacarpal_2_right", "metacarpal_3_right",
  "metacarpal_4_right", "metacarpal_5_right", "metacarpal_left", "metacarpal_left_indet", "metacarpal_1_left",
  "metacarpal_2_left", "metacarpal_3_left", "metacarpal_4_left", "metacarpal_5_left", 
  "hand_phalanx",
  "hand_phalanx_indet", 
  "hand_phalanx_proximal", "hand_phalanx_proximal_indet", "hand_phalanx_right",
  "hand_phalanx_proximal_right_indet", "hand_phalanx_proximal_1_right", "hand_phalanx_proximal_2_right",
  "hand_phalanx_proximal_3_right", "hand_phalanx_proximal_4_right", "hand_phalanx_proximal_5_right",
  "hand_phalanx_left", "hand_phalanx_proximal_left_indet", "hand_phalanx_proximal_1_left",
  "hand_phalanx_proximal_2_left", "hand_phalanx_proximal_3_left", "hand_phalanx_proximal_4_left",
  "hand_phalanx_proximal_5_left", 
  "hand_phalanx_middle", "hand_phalanx_middle_indet", "hand_middle_right",
  "hand_phalanx_middle_right_indet", "hand_phalanx_middle_2_right", "hand_phalanx_middle_3_right",
  "hand_phalanx_middle_4_right", "hand_phalanx_middle_5_right", "hand_middle_left", "hand_phalanx_middle_left_indet",
  "hand_phalanx_middle_2_left", "hand_phalanx_middle_3_left", "hand_phalanx_middle_4_left",
  "hand_phalanx_middle_5_left", 
  "hand_phalanx_distal", "hand_phalanx_distal_indet", "hand_phalanx_distal_right",
  "hand_phalanx_distal_right_indet", "hand_phalanx_distal_right_1", "hand_phalanx_distal_right_2",
  "hand_phalanx_distal_right_3", "hand_phalanx_distal_right_4", "hand_phalanx_distal_right_5", "hand_distal_left",
  "hand_phalanx_distal_left_indet", "hand_phalanx_distal_left_1", "hand_phalanx_distal_left_2",
  "hand_phalanx_distal_left_3", "hand_phalanx_distal_left_4", "hand_phalanx_distal_left_5", 
  
  "pelvic", "pelvic_indet",
  "sacrum", "sacrum_indet", "s1", "s2", "s3", "s4", "s5", "coccyx", "coccyx_indet", "cx1", "cx2", "cx3", "cx4",
  "os_coxae", "os_coxae_indet", "os_coxae_right", "os_coxae_right_indet", "ilium_right", "ischium_right",
  "pubis_right", "os_coxae_left", "os_coxae_left_indet", "ilium_left", "ischium_left", "pubis_left", 
  
  "legs",
  "legs_indet", "femur", "femur_indet", "femur_right", "femur_left", "tibia", "tibia_indet", "tibia_right",
  "tibia_left", "fibula", "fibula_indet", "fibula_right", "fibula_left", "patella", "patella_indet", "patella_right",
  "patella_left", 
  
  "foot", "foot_indet", 
  "tarsal", "tarsal_indet", "tarsal_right", "tarsal_right_indet", "talus_right",
  "calcaneous_right", "cuboid_right", "navicular_right", "lateral_cuneiform_right", "intermediate_cuneiform_right",
  "medial_cuneiform_right", "tarsal_left", "tarsal_left_indet", "talus_left", "calcaneous_left", "cuboid_left",
  "navicular_left", "lateral_cuneiform_left", "intermediate_cuneiform_left", "medial_cuneiform_left", 
  "metatarsal",
  "metatarsal_indet", "metatarsal_right", "metatarsal_right_indet", "metatarsal_1_right", "metatarsal_2_right",
  "metatarsal_3_right", "metatarsal_4_right", "metatarsal_5_right", "metatarsal_left", "metatarsal_left_indet",
  "metatarsal_1_left", "metatarsal_2_left", "metatarsal_3_left", "metatarsal_4_left", "metatarsal_5_left",
  "foot_phalanx", "foot_phalanx_indet", 
  "foot_phalanx_proximal", "foot_phalanx_proximal_indet", "foot_proximal_right",
  "foot_phalanx_proximal_right_indet", "foot_phalanx_proximal_1_right", "foot_phalanx_proximal_2_right",
  "foot_phalanx_proximal_3_right", "foot_phalanx_proximal_4_right", "foot_phalanx_proximal_5_right",
  "foot_phalanx_left", "foot_phalanx_proximal_left_indet", "foot_phalanx_proximal_1_left",
  "foot_phalanx_proximal_2_left", "foot_phalanx_proximal_3_left", "foot_phalanx_proximal_4_left",
  "foot_phalanx_proximal_5_left", 
  "foot_phalanx_middle", "foot_phalanx_middle_indet", "foot_middle_right",
  "foot_phalanx_middle_right_indet", "foot_phalanx_middle_2_right", "foot_phalanx_middle_3_right",
  "foot_phalanx_middle_4_right", "foot_phalanx_middle_5_right", "foot_middle_left", "foot_phalanx_middle_left_indet",
  "foot_phalanx_middle_2_left", "foot_phalanx_middle_3_left", "foot_phalanx_middle_4_left",
  "foot_phalanx_middle_5_left", 
  "foot_phalanx_distal", "foot_phalanx_distal_indet", "foot_phalanx_distal_right",
  "foot_phalanx_distal_right_indet", "foot_phalanx_distal_right_1", "foot_phalanx_distal_right_2",
  "foot_phalanx_distal_right_3", "foot_phalanx_distal_right_4", "foot_phalanx_distal_right_5", "foot_distal_left",
  "foot_phalanx_distal_left_indet", "foot_phalanx_distal_left_1", "foot_phalanx_distal_left_2",
  "foot_phalanx_distal_left_3", "foot_phalanx_distal_left_4", "foot_phalanx_distal_left_5"
)
Level <- c(
  #cranium
  1, 2, 3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 3, 3, 3, 
  2, 3, 3, 3, 3, 4, 4, 4, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5, 4, 5, 5, 5, 4, 5,
  5, 5, 
  3, 
  
  #teeth
  1, 2, 3, #3, 
  4, 4, 5, 5, 5, 4, 5, 5, 5, 3, 4, 4, 5, 5, 5, 4, 5, 5, 5, 3, 4, 4, 5, 5, 5, 4, 5, 5, 5, 3, 4, 4, 5,
  5, 5, 4, 5, 5, 5, 3, 4, 4, 5, 5, 5, 4, 5, 5, 5, 
  2, 3, #3, 
  4, 4, 5, 5, 5, 4, 5, 5, 5, 3, 4, 4, 5, 5, 5, 4, 5, 5, 5, 3, 4,
  4, 5, 5, 5, 4, 5, 5, 5, 3, 4, 4, 5, 5, 5, 4, 5, 5, 5, 3, 4, 4, 5, 5, 5, 4, 5, 5, 5, 3, 4, 4, 5, 5, 5, 4, 5, 5, 5, 3,
  4, 4, 5, 5, 5, 4, 5, 5, 5, 3, 4, 4, 5, 5, 5, 4, 5, 5, 5, 
  
  #hyoid
  1, 
  
  #vertebrae
  1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 
  
  #thorax
  1, 2, 2, 2, 3, 3, 3, 3, 
  3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 
  3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4,  
  
  #shoulder
  1, 2, 2, 3, 3, 3, 2, 3, 3, 3, 
  
  #arms
  1, 2, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 
  
  #hands
  1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4,
  4, 4, 4, 4, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 3, 4, 4, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 3, 4, 4,
  5, 5,  5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 
  
  #pelvic
  1, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 
  2, 3, 3, 4, 4, 4, 4, 3, 4, 4, 4, 4, 
  
  #legs
  1, 2, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 
  
  #foot
  1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4,
  4, 4, 2, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 3, 4, 4,
  5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 3, 4, 4, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5
)

max_Level <- c(
  #cranium
  1, 2, 5, 5, 3, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 3, 5, 4, 5, 5, 5, 4, 5, 5, 5, 4, 5, 5,
  5, 5, 
  
  #teeth
  1, 2, #5, 
  3, 5, 4, 5, 5, 5, 4, 5, 5, 5, 3, 5, 4, 5, 5, 5, 4, 5, 5, 5, 3, 5, 4, 5, 5, 5, 4, 5, 5, 5, 3, 5, 4, 5, 5,
  5, 4, 5, 5, 5, 3, 5, 4, 5, 5, 5, 4, 5, 5, 5, 
  2, #5, 
  3, 5, 4, 5, 5, 5, 4, 5, 5, 5, 3, 5, 4, 5, 5, 5, 4, 5, 5, 5, 3, 5, 4,
  5, 5, 5, 4, 5, 5, 5, 3, 5, 4, 5, 5, 5, 4, 5, 5, 5, 3, 5, 4, 5, 5, 5, 4, 5, 5, 5, 3, 5, 4, 5, 5, 5, 4, 5, 5, 5, 3, 5, 4,
  5, 5, 5, 4, 5, 5, 5, 3, 5, 4, 5, 5, 5, 4, 5, 5, 5, 
  
  #hyoid
  1, 
  
 #vertebrae
  1, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 
  
  #thorax
  1, 4, 4, 2, 4, 4, 4, 4, 
  3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 
  3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 
  
  #shoulder
  1, 3, 2, 3, 3, 3, 2, 3, 3, 3, 
  
  #arms
  1, 3, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 
  
  #hands
  1, 5, 2, 5,
  3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 5, 3, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 2, 5, 3, 5,
  4, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 3, 5, 4, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 3, 5, 4, 5, 5, 5, 5, 5, 5, 4, 5, 5,
  5, 5, 5, 5, 
  
  #pelvic
  1, 4, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 
  2, 4, 3, 4, 4, 4, 4, 3, 4, 4, 4, 4, 
  
  #legs
  1, 3, 2, 3, 3, 3, 2, 3, 3,
  3, 2, 3, 3, 3, 2, 3, 3, 3, 
  
  #foot
  1, 5, 2, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 2, 5, 3, 5, 5, 5, 5, 5, 5,
  3, 5, 5, 5, 5, 5, 5, 2, 5, 3, 5, 4, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 3, 5, 4, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 3,
  5, 4, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5
)


Bone_Count <- c(
  #cranium
  1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 2, 
  6, 6, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 
  #teeth
  1, 20, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 
  4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 
  32, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 
  4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 
  4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 
  #hyoid
  1, 
  #vertebrae
  24, 24, 7, 7, 1, 1, 1, 1, 1, 1, 1, 12, 12, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1, 
  #thorax
  25, 25, 1, 24, 72, 24, 12, 12, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 
  2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 
  #shoulder
  4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 
  #arms
  6, 6, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 
  #hands
  54, 54, 16, 16, 8, 8, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 1, 1, 1, 1, 1, 1, 1, 1, 
  10, 10, 5, 5, 1, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1, 
  8, 28, 10, 10, 5, 5, 1, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1, 8, 8, 4, 4, 1, 1, 1, 1, 4, 4, 1, 1, 1, 1, 
  10, 10, 5, 5, 1, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1, 
   #pelvic
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    #legs
  8, 8, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 
    #foot
  52, 52, 14, 14, 7, 7, 1, 1, 1, 1, 1, 1, 1, 7, 7, 1, 1, 1, 1, 1, 1, 1, 
  10, 10, 5, 5, 1, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1, 
  28, 28, 10, 10, 5, 5, 1, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1, 8, 8, 4, 4, 1, 1, 1, 1, 4, 4, 1, 1, 1, 1, 
  10, 10, 5, 5, 1, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1
)
Skeleton_Section <- c(
  "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium",
  "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium",
  "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium",
  "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium", "cranium",
  "cranium", "cranium", 
  
  #"teeth", "teeth",
  "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth",
  "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth",
  "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth",
  "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth",
  "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth",
  "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth",
  "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth",
  "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth",
  "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth",
  "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", "teeth",
  "teeth", "teeth", "teeth", "teeth", "teeth", "teeth", 
  
  "hyoid", 
  
  "vertebrae", "vertebrae", "vertebrae", "vertebrae",
  "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae",
  "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae",
  "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae", "vertebrae",
  "vertebrae", 
  
  "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", 
  "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax",
  "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax",
  "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax", "thorax",
  "thorax", 
  
  "shoulder", "shoulder", "shoulder", "shoulder", "shoulder", "shoulder", "shoulder", "shoulder", "shoulder",
  "shoulder", 
  
  "arms", "arms", "arms", "arms", "arms", "arms", "arms", "arms", "arms", "arms", "arms", "arms", "arms",
  "arms", 
  
  "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands",
  "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands",
  "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands",
  "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands",
  "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands",
  "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands",
  "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", "hands", 
  
  "pelvic",
  "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic",
  "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic", "pelvic",
  "pelvic", "pelvic", "pelvic", "pelvic", 
  
  "legs", "legs", "legs", "legs", "legs", "legs", "legs", "legs", "legs",
  "legs", "legs", "legs", "legs", "legs", "legs", "legs", "legs", "legs", 
  
  "foot", "foot", "foot", "foot", "foot",
  "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot",
  "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot",
  "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot",
  "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot",
  "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot",
  "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot", "foot"
)

elements_data <- data.frame(
  ElementsOne = ElementsOne,
  Level = Level,
  max_Level = max_Level,
  Bone_Count = Bone_Count,
  Skeleton_Section = Skeleton_Section,
  stringsAsFactors = FALSE
)


# Ensure types
elements_data$Level <- as.integer(elements_data$Level)
elements_data$Bone_Count <- as.integer(elements_data$Bone_Count)
elements_data$ElementsOne <- as.character(elements_data$ElementsOne)
elements_data$Skeleton_Section <- as.character(elements_data$Skeleton_Section)

# Section order (xlsx includes teeth)
sections_order <- c("cranium","teeth","hyoid","vertebrae","thorax","shoulder","arms","hands","pelvic","legs","foot")
elements_data$Skeleton_Section <- factor(elements_data$Skeleton_Section, levels = sections_order)



####### mappings ####
category_mappings = list(



#### CRANIUM 
# Level 2 -> Level 1
category_mappings_cranium_2_to_1 <- list(
  cranium = c("neurocranium_bone", "splanchnocranium")
),

# Level 3 -> Level 1
category_mappings_cranium_3_to_1 <- list(
  cranium = c(
    # neurocranium (L3)
    "neurocranium_bone_indet",
    "frontal_bone", "parietal_bone", "temporal", "occipital", "sphenoid", "ethmoid",

    # splanchnocranium (L3)
    "splanchnocranium_indet",
    "maxilla", "palatine_bone", "zygomatic_bone", "nasal_bone", "lacrimal_bone",
    "vomer_bone", "inferior_nasal_concha",

    # other L3
    "ossicles_of_the_middle_ear",
    "mandible"
  )
),

# Level 4 -> Level 1
category_mappings_cranium_4_to_1 <- list(
  cranium = c(
    "neurocranium_bone_indet","frontal_bone",
    "parietal_bone_indet", "parietal_bone_right", "parietal_bone_left",
    "temporal_bone_indet", "temporal_right", "temporal_left",
    "occipital", "sphenoid", "ethmoid",
    "splanchnocranium_indet",    
    "maxilla", "palatine_bone",
    "zygomatic_bone_indet", "zygomatic_bone_right", "zygomatic_bone_left",
    "nasal_bone", "lacrimal_bone", "vomer_bone", "inferior_nasal_concha",

    "ossicles_of_the_middle_ear_indet",
    "malleus", "incus", "stapes",

    "mandible"
  )
),

# Level 5 -> Level 1
category_mappings_cranium_5_to_1 <- list(
  cranium = c("neurocranium_bone_indet",
    "frontal_bone",
    "parietal_bone_indet", "parietal_bone_right", "parietal_bone_left",
    "temporal_bone_indet", "temporal_right", "temporal_left",
    "occipital", "sphenoid", "ethmoid",
    "splanchnocranium_indet",
    "maxilla", "palatine_bone",
    "zygomatic_bone_indet", "zygomatic_bone_right", "zygomatic_bone_left",
    "nasal_bone", "lacrimal_bone", "vomer_bone", "inferior_nasal_concha",

    "ossicles_of_the_middle_ear_indet",
    "malleus_indet", "malleus_right", "malleus_left",
    "incus_indet", "incus_right", "incus_left",
    "stapes_indet", "stapes_right", "stapes_left",

    "mandible"
  )
),

# Level 3 -> Level 2  (ONLY things that actually belong to L2 groups)
category_mappings_cranium_3_to_2 <- list(
  neurocranium_bone = c(
    "neurocranium_bone_indet",
    "frontal_bone", "parietal_bone", "temporal", "occipital", "sphenoid", "ethmoid"
  ),
  splanchnocranium = c(
    "splanchnocranium_indet",
    "maxilla", "palatine_bone", "zygomatic_bone", "nasal_bone", "lacrimal_bone",
    "vomer_bone", "inferior_nasal_concha",
    "ossicles_of_the_middle_ear","mandible"

  )
),

# Level 4 -> Level 2
category_mappings_cranium_4_to_2 <- list(
  neurocranium_bone = c("neurocranium_bone_indet",
    "frontal_bone",
    "parietal_bone_indet", "parietal_bone_right", "parietal_bone_left",
    "temporal_bone_indet", "temporal_right", "temporal_left",
    "occipital", "sphenoid", "ethmoid"
  ),
  splanchnocranium = c("splanchnocranium_indet",
    "maxilla", "palatine_bone",
    "zygomatic_bone_indet", "zygomatic_bone_right", "zygomatic_bone_left",
    "nasal_bone", "lacrimal_bone", "vomer_bone", "inferior_nasal_concha",
    "ossicles_of_the_middle_ear_indet","malleus","incus","stapes","mandible"
  )
),

# Level 5 -> Level 2
category_mappings_cranium_5_to_2 <- list(
    neurocranium_bone = c("neurocranium_bone_indet",
      "frontal_bone",
      "parietal_bone_indet", "parietal_bone_right", "parietal_bone_left",
      "temporal_bone_indet", "temporal_right", "temporal_left",
      "occipital", "sphenoid", "ethmoid"
    ),
    splanchnocranium = c("splanchnocranium_indet",
      "maxilla", "palatine_bone",
      "zygomatic_bone_indet", "zygomatic_bone_right", "zygomatic_bone_left",
      "nasal_bone", "lacrimal_bone", "vomer_bone", "inferior_nasal_concha",
      "ossicles_of_the_middle_ear_indet","malleus_indet","malleus_right",
      "malleus_left","incus_indet","incus_right","incus_left","stapes_indet",
      "stapes_right","stapes_left",      
      "mandible"
    )
  ),

# Level 4 -> Level 3
category_mappings_cranium_4_to_3 <- list(
    neurocranium_bone_indet  = c("neurocranium_bone_indet"),  
  frontal_bone  = c("frontal_bone"),
  parietal_bone = c("parietal_bone_indet", "parietal_bone_right", "parietal_bone_left"),
  temporal      = c("temporal_bone_indet", "temporal_right", "temporal_left"),
  occipital     = c("occipital"),
  sphenoid      = c("sphenoid"),
  ethmoid       = c("ethmoid"),

  splanchnocranium_indet  = c("splanchnocranium_indet"), 
  maxilla       = c("maxilla"),
  palatine_bone = c("palatine_bone"),
  zygomatic_bone = c("zygomatic_bone_indet", "zygomatic_bone_right", "zygomatic_bone_left"),
  nasal_bone    = c("nasal_bone"),
  lacrimal_bone = c("lacrimal_bone"),
  vomer_bone    = c("vomer_bone"),
  inferior_nasal_concha = c("inferior_nasal_concha"),

  ossicles_of_the_middle_ear = c("ossicles_of_the_middle_ear_indet", "malleus", "incus", "stapes"),
  mandible = c("mandible")
),

# Level 5 -> Level 3
category_mappings_cranium_5_to_3 <- list(
    neurocranium_bone_indet  = c("neurocranium_bone_indet"),  
  frontal_bone  = c("frontal_bone"),
  parietal_bone = c("parietal_bone_indet", "parietal_bone_right", "parietal_bone_left"),
  temporal      = c("temporal_bone_indet", "temporal_right", "temporal_left"),
  occipital     = c("occipital"),
  sphenoid      = c("sphenoid"),
  ethmoid       = c("ethmoid"),

  splanchnocranium_indet  = c("splanchnocranium_indet"), 
  maxilla       = c("maxilla"),
  palatine_bone = c("palatine_bone"),
  zygomatic_bone = c("zygomatic_bone_indet", "zygomatic_bone_right", "zygomatic_bone_left"),
  nasal_bone    = c("nasal_bone"),
  lacrimal_bone = c("lacrimal_bone"),
  vomer_bone    = c("vomer_bone"),
  inferior_nasal_concha = c("inferior_nasal_concha"),

  ossicles_of_the_middle_ear = c(
    "ossicles_of_the_middle_ear_indet",
    "malleus_indet", "malleus_right", "malleus_left",
    "incus_indet", "incus_right", "incus_left",
    "stapes_indet", "stapes_right", "stapes_left"
  ),
  mandible = c("mandible")
),

# Level 5 -> Level 4
category_mappings_cranium_5_to_4 <- list(
    neurocranium_bone_indet  = c("neurocranium_bone_indet"),
  frontal_bone = c("frontal_bone"),
  parietal_bone_indet = c("parietal_bone_indet"),
  parietal_bone_right = c("parietal_bone_right"),
  parietal_bone_left  = c("parietal_bone_left"),
  temporal_bone_indet = c("temporal_bone_indet"),
  temporal_right      = c("temporal_right"),
  temporal_left       = c("temporal_left"),
  occipital = c("occipital"),
  sphenoid  = c("sphenoid"),
  ethmoid   = c("ethmoid"),

  splanchnocranium_indet  = c("splanchnocranium_indet"), 
  maxilla       = c("maxilla"),
  palatine_bone = c("palatine_bone"),
  zygomatic_bone_indet = c("zygomatic_bone_indet"),
  zygomatic_bone_right = c("zygomatic_bone_right"),
  zygomatic_bone_left  = c("zygomatic_bone_left"),
  nasal_bone    = c("nasal_bone"),
  lacrimal_bone = c("lacrimal_bone"),
  vomer_bone    = c("vomer_bone"),
  inferior_nasal_concha = c("inferior_nasal_concha"),

  ossicles_of_the_middle_ear_indet = c("ossicles_of_the_middle_ear_indet"),
  malleus = c("malleus_indet", "malleus_right", "malleus_left"),
  incus   = c("incus_indet", "incus_right", "incus_left"),
  stapes  = c("stapes_indet", "stapes_right", "stapes_left"),
  
  mandible = c("mandible")
),



#Hyoid 1


#### TEETH
# Level 2 -> Level 1
category_mappings_teeth_2_to_1 <- list(
  teeth = c("deciduous", "permanent")
  #teeth = c("deciduous_indet", "permanent")
),

# Level 3 -> Level 1
category_mappings_teeth_3_to_1 <- list(
  teeth = c(
      #"deciduous_indet", 
      "di1", "di2", "dc", "dm1", "dm2", 
      #"permanent_indet", 
      "i1", "i2", "c", "pm3", "pm4",
      "m1", "m2", "m3"
    )
)
,
# Level 4 -> Level 1
category_mappings_teeth_4_to_1 <- list(
  teeth = c(
      #"deciduous_indet", 
      "dl1_indet", "di1_superior", "dl1_inferior", "dl2_indet", "di2_superior",
      "dl2_inferior", "dc_indet", "dc_superior", "dc_inferior", "dm1_indet", "dm1_superior", "dm1_inferior",
      "dm2_indet", "dm2_superior", "dm2_inferior", 
      #"permanent_indet", 
      "i1_indet", "i1_superior",
      "i1_inferior", "i2_indet", "i2_superior", "i2_inferior", "c_indet", "c_superior", "c_inferior",
      "pm3_indet", "pm3_superior", "pm3_inferior", "pm4_indet", "pm4_superior", "pm4_inferior", "m1_indet",
      "m1_superior", "m1_inferior", "m2_indet", "m2_superior", "m2_inferior", "m3_indet", "m3_superior",
      "m3_inferior"
    )
)
,
# Level 5 -> Level 1
category_mappings_teeth_5_to_1 <- list(
  teeth = c(
      #"deciduous_indet", 
      "dl1_indet", "di1_superior_indet", "di1_superior_right", "di1_superior_left",
      "dl1_inferior_indet", "dl1_inferior_right", "dl1_inferior_left", "dl2_indet", "di2_superior_indet",
      "di2_superior_right", "di2_superior_left", "di2_inferior_indet", "di2_inferior_right",
      "di2_inferior_left", "dc_indet", "dc_superior_indet", "dc_superior_right", "dc_superior_left",
      "dc_inferior_indet", "dc_inferior_right", "dc_inferior_left", "dm1_indet", "dm1_superior_indet",
      "dm1_superior_left", "dm1_superior_right", "dm1_inferior_indet", "dm1_inferior_right",
      "dm1_inferior_left", "dm2_indet", "dm2_superior_indet", "dm2_superior_right", "dm2_superior_left",
      "dm2_inferior_indet", "dm2_inferior_right", "dm2_inferior_left", 
      #"permanent_indet", 
      "i1_indet",
      "i1_superior_indet", "i1_superior_right", "i1_superior_left", "i1_inferior_indet", "i1_inferior_right",
      "i1_inferior_left", "i2_indet", "i2_superior_indet", "i2_superior_right", "i2_superior_left",
      "i2_inferior_indet", "i2_inferior_right", "i2_inferior_left", "c_indet", "c_superior_indet",
      "c_superior_right", "c_superior_left", "c_inferior_indet", "c_inferior_right", "c_inferior_left",
      "pm3_indet", "pm3_superior_indet", "pm3_superior_right", "pm3_superior_left", "pm3_inferior_indet",
      "pm3_inferior_right", "pm3_inferior_left", "pm4_indet", "pm4_superior_indet", "pm4_superior_right",
      "pm4_superior_left", "pm4_inferior_indet", "pm4_inferior_right", "pm4_inferior_left", "m1_indet",
      "m1_superior_indet", "m1_superior_right", "m1_superior_left", "m1_inferior_indet", "m1_inferior_right",
      "m1_inferior_left", "m2_indet", "m2_superior_indet", "m2_superior_right", "m2_superior_left",
      "m2_inferior_indet", "m2_inferior_right", "m2_inferior_left", "m3_indet", "m3_superior_indet",
      "m3_superior_right", "m3_superior_left", "m3_inferior_indet", "m3_inferior_right", "m3_inferior_left"
    )
),

# Level 3 -> Level 2
category_mappings_teeth_3_to_2 <- list(
  #deciduous = c("deciduous_indet", "di1", "di2", "dc", "dm1", "dm2"),
  deciduous = c("di1", "di2", "dc", "dm1", "dm2"),
  permanent = c( #"permanent_indet", 
                "i1", "i2", "c", "pm3", "pm4", "m1", "m2", "m3")
),

# Level 4 -> Level 2
category_mappings_teeth_4_to_2 <- list(
  deciduous = c( 
        #"deciduous_indet", 
        "dl1_indet", "di1_superior", "dl1_inferior", "dl2_indet", "di2_superior",
        "dl2_inferior", "dc_indet", "dc_superior", "dc_inferior", "dm1_indet", "dm1_superior", "dm1_inferior",
        "dm2_indet", "dm2_superior", "dm2_inferior"
      ),
  permanent = c(
        #"permanent_indet", 
        "i1_indet", "i1_superior", "i1_inferior", "i2_indet", "i2_superior", "i2_inferior",
        "c_indet", "c_superior", "c_inferior", "pm3_indet", "pm3_superior", "pm3_inferior", "pm4_indet",
        "pm4_superior", "pm4_inferior", "m1_indet", "m1_superior", "m1_inferior", "m2_indet", "m2_superior",
        "m2_inferior", "m3_indet", "m3_superior", "m3_inferior"
      )
),

# Level 5 -> Level 2
category_mappings_teeth_5_to_2 <- list(
  deciduous = c(
        #"deciduous_indet", 
        "dl1_indet", "di1_superior_indet", "di1_superior_right", "di1_superior_left",
        "dl1_inferior_indet", "dl1_inferior_right", "dl1_inferior_left", "dl2_indet", "di2_superior_indet",
        "di2_superior_right", "di2_superior_left", "di2_inferior_indet", "di2_inferior_right",
        "di2_inferior_left", "dc_indet", "dc_superior_indet", "dc_superior_right", "dc_superior_left",
        "dc_inferior_indet", "dc_inferior_right", "dc_inferior_left", "dm1_indet", "dm1_superior_indet",
        "dm1_superior_left", "dm1_superior_right", "dm1_inferior_indet", "dm1_inferior_right",
        "dm1_inferior_left", "dm2_indet", "dm2_superior_indet", "dm2_superior_right", "dm2_superior_left",
        "dm2_inferior_indet", "dm2_inferior_right", "dm2_inferior_left"
      ),
  permanent = c(
        #"permanent_indet", 
        "i1_indet", "i1_superior_indet", "i1_superior_right", "i1_superior_left",
        "i1_inferior_indet", "i1_inferior_right", "i1_inferior_left", "i2_indet", "i2_superior_indet",
        "i2_superior_right", "i2_superior_left", "i2_inferior_indet", "i2_inferior_right", "i2_inferior_left",
        "c_indet", "c_superior_indet", "c_superior_right", "c_superior_left", "c_inferior_indet",
        "c_inferior_right", "c_inferior_left", "pm3_indet", "pm3_superior_indet", "pm3_superior_right",
        "pm3_superior_left", "pm3_inferior_indet", "pm3_inferior_right", "pm3_inferior_left", "pm4_indet",
        "pm4_superior_indet", "pm4_superior_right", "pm4_superior_left", "pm4_inferior_indet",
        "pm4_inferior_right", "pm4_inferior_left", "m1_indet", "m1_superior_indet", "m1_superior_right",
        "m1_superior_left", "m1_inferior_indet", "m1_inferior_right", "m1_inferior_left", "m2_indet",
        "m2_superior_indet", "m2_superior_right", "m2_superior_left", "m2_inferior_indet",
        "m2_inferior_right", "m2_inferior_left", "m3_indet", "m3_superior_indet", "m3_superior_right",
        "m3_superior_left", "m3_inferior_indet", "m3_inferior_right", "m3_inferior_left"
      )
),

# Level 4 -> Level 3
category_mappings_teeth_4_to_3 <- list(
  #deciduous_indet = c("deciduous_indet"),
  di1 = c("dl1_indet", "di1_superior", "dl1_inferior"),
  di2 = c("dl2_indet", "di2_superior", "dl2_inferior"),
  dc = c("dc_indet", "dc_superior", "dc_inferior"),
  dm1 = c("dm1_indet", "dm1_superior", "dm1_inferior"),
  dm2 = c("dm2_indet", "dm2_superior", "dm2_inferior"),
  #permanent_indet = c("permanent_indet"),
  i1 = c("i1_indet", "i1_superior", "i1_inferior"),
  i2 = c("i2_indet", "i2_superior", "i2_inferior"),
  c = c("c_indet", "c_superior", "c_inferior"),
  pm3 = c("pm3_indet", "pm3_superior", "pm3_inferior"),
  pm4 = c("pm4_indet", "pm4_superior", "pm4_inferior"),
  m1 = c("m1_indet", "m1_superior", "m1_inferior"),
  m2 = c("m2_indet", "m2_superior", "m2_inferior"),
  m3 = c("m3_indet", "m3_superior", "m3_inferior")
),

# Level 5 -> Level 3
category_mappings_teeth_5_to_3 <- list(
  #deciduous_indet = c("deciduous_indet"),
  di1 = c(
        "dl1_indet", "di1_superior_indet", "di1_superior_right", "di1_superior_left", "dl1_inferior_indet",
        "dl1_inferior_right", "dl1_inferior_left"
      ),
  di2 = c(
        "dl2_indet", "di2_superior_indet", "di2_superior_right", "di2_superior_left", "di2_inferior_indet",
        "di2_inferior_right", "di2_inferior_left"
      ),
  dc = c(
        "dc_indet", "dc_superior_indet", "dc_superior_right", "dc_superior_left", "dc_inferior_indet",
        "dc_inferior_right", "dc_inferior_left"
      ),
  dm1 = c(
        "dm1_indet", "dm1_superior_indet", "dm1_superior_left", "dm1_superior_right", "dm1_inferior_indet",
        "dm1_inferior_right", "dm1_inferior_left"
      ),
  dm2 = c(
        "dm2_indet", "dm2_superior_indet", "dm2_superior_right", "dm2_superior_left", "dm2_inferior_indet",
        "dm2_inferior_right", "dm2_inferior_left"
      ),
  #permanent_indet = c("permanent_indet"),
  i1 = c(
        "i1_indet", "i1_superior_indet", "i1_superior_right", "i1_superior_left", "i1_inferior_indet",
        "i1_inferior_right", "i1_inferior_left"
      ),
  i2 = c(
        "i2_indet", "i2_superior_indet", "i2_superior_right", "i2_superior_left", "i2_inferior_indet",
        "i2_inferior_right", "i2_inferior_left"
      ),
  c = c(
        "c_indet", "c_superior_indet", "c_superior_right", "c_superior_left", "c_inferior_indet",
        "c_inferior_right", "c_inferior_left"
      ),
  pm3 = c(
        "pm3_indet", "pm3_superior_indet", "pm3_superior_right", "pm3_superior_left", "pm3_inferior_indet",
        "pm3_inferior_right", "pm3_inferior_left"
      ),
  pm4 = c(
        "pm4_indet", "pm4_superior_indet", "pm4_superior_right", "pm4_superior_left", "pm4_inferior_indet",
        "pm4_inferior_right", "pm4_inferior_left"
      ),
  m1 = c(
        "m1_indet", "m1_superior_indet", "m1_superior_right", "m1_superior_left", "m1_inferior_indet",
        "m1_inferior_right", "m1_inferior_left"
      ),
  m2 = c(
        "m2_indet", "m2_superior_indet", "m2_superior_right", "m2_superior_left", "m2_inferior_indet",
        "m2_inferior_right", "m2_inferior_left"
      ),
  m3 = c(
        "m3_indet", "m3_superior_indet", "m3_superior_right", "m3_superior_left", "m3_inferior_indet",
        "m3_inferior_right", "m3_inferior_left"
      )
),

# Level 5 -> Level 4
category_mappings_teeth_5_to_4 <- list(
  #deciduous_indet = c("deciduous_indet"),
  dl1_indet = c("dl1_indet"),
  di1_superior = c("di1_superior_indet", "di1_superior_right", "di1_superior_left"),
  dl1_inferior = c("dl1_inferior_indet", "dl1_inferior_right", "dl1_inferior_left"),
  dl2_indet = c("dl2_indet"),
  di2_superior = c("di2_superior_indet", "di2_superior_right", "di2_superior_left"),
  dl2_inferior = c("di2_inferior_indet", "di2_inferior_right", "di2_inferior_left"),
  dc_indet = c("dc_indet"),
  dc_superior = c("dc_superior_indet", "dc_superior_right", "dc_superior_left"),
  dc_inferior = c("dc_inferior_indet", "dc_inferior_right", "dc_inferior_left"),
  dm1_indet = c("dm1_indet"),
  dm1_superior = c("dm1_superior_indet", "dm1_superior_left", "dm1_superior_right"),
  dm1_inferior = c("dm1_inferior_indet", "dm1_inferior_right", "dm1_inferior_left"),
  dm2_indet = c("dm2_indet"),
  dm2_superior = c("dm2_superior_indet", "dm2_superior_right", "dm2_superior_left"),
  dm2_inferior = c("dm2_inferior_indet", "dm2_inferior_right", "dm2_inferior_left"),
  #permanent_indet = c("permanent_indet"),
  i1_indet = c("i1_indet"),
  i1_superior = c("i1_superior_indet", "i1_superior_right", "i1_superior_left"),
  i1_inferior = c("i1_inferior_indet", "i1_inferior_right", "i1_inferior_left"),
  i2_indet = c("i2_indet"),
  i2_superior = c("i2_superior_indet", "i2_superior_right", "i2_superior_left"),
  i2_inferior = c("i2_inferior_indet", "i2_inferior_right", "i2_inferior_left"),
  c_indet = c("c_indet"),
  c_superior = c("c_superior_indet", "c_superior_right", "c_superior_left"),
  c_inferior = c("c_inferior_indet", "c_inferior_right", "c_inferior_left"),
  pm3_indet = c("pm3_indet"),
  pm3_superior = c("pm3_superior_indet", "pm3_superior_right", "pm3_superior_left"),
  pm3_inferior = c("pm3_inferior_indet", "pm3_inferior_right", "pm3_inferior_left"),
  pm4_indet = c("pm4_indet"),
  pm4_superior = c("pm4_superior_indet", "pm4_superior_right", "pm4_superior_left"),
  pm4_inferior = c("pm4_inferior_indet", "pm4_inferior_right", "pm4_inferior_left"),
  m1_indet = c("m1_indet"),
  m1_superior = c("m1_superior_indet", "m1_superior_right", "m1_superior_left"),
  m1_inferior = c("m1_inferior_indet", "m1_inferior_right", "m1_inferior_left"),
  m2_indet = c("m2_indet"),
  m2_superior = c("m2_superior_indet", "m2_superior_right", "m2_superior_left"),
  m2_inferior = c("m2_inferior_indet", "m2_inferior_right", "m2_inferior_left"),
  m3_indet = c("m3_indet"),
  m3_superior = c("m3_superior_indet", "m3_superior_right", "m3_superior_left"),
  m3_inferior = c("m3_inferior_indet", "m3_inferior_right", "m3_inferior_left")
),



#### HYOID


#### VERTEBRAE

# Level 2 -> Level 1
category_mappings_vertebrae_2_to_1 <- list(
  vertebrae = c("vertebrae_indet", "cervical", "thoracic", "lumbar")
),

# Level 3 -> Level 1
category_mappings_vertebrae_3_to_1 <- list(
  vertebrae = c(
      "vertebrae_indet", "cervical_indet", "c1_atlas", "c2_axis", "c3", "c4", "c5", "c6",
      "c7", "thoracic_indet", "t1", "t2", "t3", "t4", "t5", "t6",
      "t7", "t8", "t9", "t10", "t11", "t12", "lumbar_indet", "l1",
      "l2", "l3", "l4", "l5"
    )
),

# Level 3 -> Level 2
category_mappings_vertebrae_3_to_2 <- list(
  vertebrae_indet = c("vertebrae_indet"),
  cervical = c("cervical_indet", "c1_atlas", "c2_axis", "c3", "c4", "c5", "c6", "c7"),
  thoracic = c(
      "thoracic_indet", "t1", "t2", "t3", "t4", "t5", "t6", "t7",
      "t8", "t9", "t10", "t11", "t12"
    ),
  lumbar = c("lumbar_indet", "l1", "l2", "l3", "l4", "l5")
),


#### THORAX

# Level 2 -> Level 1
category_mappings_thorax_2_to_1 <- list(
  thorax = c("thorax_indet", "sternum", "ribs" )
),

# Level 3 -> Level 1
category_mappings_thorax_3_to_1 <- list(
  thorax = c( "thorax_indet", "sternum", "rib_indet", "rib_fragments", "ribs_right", "ribs_left",  
  "rib_1",  "rib_2", "rib_3", "rib_4", "rib_5",  "rib_6", "rib_7",  
  "rib_8", "rib_9", "rib_10", "rib_11", "rib_12"
)
),

# Level 4 -> Level 1
category_mappings_thorax_4_to_1 <- list(
  thorax = c(
    "thorax_indet", "sternum", "rib_indet", "rib_fragments", "ribs_right", "ribs_left", 
     "rib_1_right", "rib_1_left",
     "rib_2_right", "rib_2_left",  "rib_3_right", "rib_3_left", "rib_4_right",
     "rib_4_left",  "rib_5_right", "rib_5_left",  
     "rib_6_right", "rib_6_left", 
     "rib_7_right", "rib_7_left", "rib_8_right", "rib_8_left",  "rib_9_right", "rib_9_left",
     "rib_10_right", "rib_10_left",  "rib_11_right", "rib_11_left", "rib_12_right", "rib_12_left"
  )
),

# Level 3 -> Level 2
category_mappings_thorax_3_to_2 <- list(
  thorax_indet = c("thorax_indet"),
  sternum = c("sternum"),
  ribs = c("rib_indet", "rib_fragments", "ribs_right", "ribs_left",  
  "rib_1",  "rib_2", "rib_3", "rib_4", "rib_5",  "rib_6", "rib_7",  
  "rib_8", "rib_9", "rib_10", "rib_11", "rib_12")
),

# Level 4 -> Level 2
category_mappings_thorax_4_to_2 <- list(
  thorax_indet = c("thorax_indet"),
  sternum = c("sternum"),
  ribs = c(
      "rib_indet", "rib_fragments", 
    "ribs_right", "ribs_left", 
     "rib_1_right", "rib_1_left",
     "rib_2_right", "rib_2_left",  "rib_3_right", "rib_3_left", "rib_4_right",
     "rib_4_left",  "rib_5_right", "rib_5_left",  
     "rib_6_right", "rib_6_left", 
     "rib_7_right", "rib_7_left", "rib_8_right", "rib_8_left",  "rib_9_right", "rib_9_left",
     "rib_10_right", "rib_10_left",  "rib_11_right", "rib_11_left", "rib_12_right", "rib_12_left"
    )
),

# Level 4 -> Level 3
category_mappings_thorax_4_to_3 <- list(
  thorax_indet = c("thorax_indet"),
  sternum = c("sternum"),
  rib_indet = c("rib_indet"),
  rib_fragments = c("rib_fragments"),
  ribs_right = c("ribs_right"),
  ribs_left = c( "ribs_left" ),
  rib_1 = c("rib_1_right", "rib_1_left"),
  rib_2 = c("rib_2_right", "rib_2_left"),
  rib_3 = c("rib_3_right", "rib_3_left"),
  rib_4 = c("rib_4_right", "rib_4_left"),
  rib_5 = c("rib_5_right", "rib_5_left"),
  rib_6 = c("rib_6_right", "rib_6_left"),
  rib_7 = c("rib_7_right", "rib_7_left"),
  rib_8 = c("rib_8_right", "rib_8_left"),
  rib_9 = c("rib_9_right", "rib_9_left"),
  rib_10 = c("rib_10_right", "rib_10_left"),
  rib_11 = c("rib_11_right", "rib_11_left"),
  rib_12 = c("rib_12_right", "rib_12_left")  
),


#### SHOULDER

# Level 2 -> Level 1
category_mappings_shoulder_2_to_1 <- list(
  shoulder = c("shoulder_indet", "clavicle", "scapula")
),

# Level 3 -> Level 1
category_mappings_shoulder_3_to_1 <- list(
  shoulder = c(
      "shoulder_indet", "clavicle_indet", "clavicle_right", "clavicle_left", "scapula_indet", "scapula_right", "scapula_left"
    )
),

# Level 3 -> Level 2
category_mappings_shoulder_3_to_2 <- list(
  shoulder_indet = c("shoulder_indet"),
  clavicle = c("clavicle_indet", "clavicle_right", "clavicle_left"),
  scapula = c("scapula_indet", "scapula_right", "scapula_left")
),


#### ARMS

# Level 2 -> Level 1
category_mappings_arms_2_to_1 <- list(
  arms = c("arms_indet", "humerus", "radius", "ulna")
),

# Level 3 -> Level 1
category_mappings_arms_3_to_1 <- list(
  arms = c(
      "arms_indet", "humerus_indet", "humerus_right", "humerus_left", "radius_indet", "radius_right", "radius_left", "ulna_indet",
      "ulna_right", "ulna_left"
    )
),

# Level 3 -> Level 2
category_mappings_arms_3_to_2 <- list(
  arms_indet = c("arms_indet"),
  humerus = c("humerus_indet", "humerus_right", "humerus_left"),
  radius = c("radius_indet", "radius_right", "radius_left"),
  ulna = c("ulna_indet", "ulna_right", "ulna_left")
),


#### HANDS

# Level 2 -> Level 1
category_mappings_hands_2_to_1 <- list(
  hands = c("carpal", "hand_phalanx", "hands_indet", "metacarpal")
),

# Level 3 -> Level 1
category_mappings_hands_3_to_1 <- list(
  hands = c(
    "hands_indet",
    "carpal_indet", "carpal_left", "carpal_right",
    "hand_phalanx_distal", "hand_phalanx_indet", "hand_phalanx_middle", "hand_phalanx_proximal",
    
    "metacarpal_indet", "metacarpal_left", "metacarpal_right"
  )
),

# Level 4 -> Level 1
category_mappings_hands_4_to_1 <- list(
  hands = c(
    "hands_indet",
    "carpal_indet", "carpal_left_indet", "carpal_right_indet",
    "pisiform_left", "pisiform_right",
    "scaphoid_left", "scaphoid_right",
    "trapezium_left", "trapezium_right",
    "trapezoid_left", "trapezoid_right",
    "triquetrum_left", "triquetrum_right",
    "lunate_left", "lunate_right",
    "capitate_left", "capitate_right",
    "hamate_left", "hamate_right",
    
    "hand_distal_left", "hand_middle_left", "hand_middle_right",
    "hand_phalanx_distal_indet", "hand_phalanx_distal_right",
    "hand_phalanx_indet", "hand_phalanx_left",
    "hand_phalanx_middle_indet", "hand_phalanx_proximal_indet",
    "hand_phalanx_right",
    
    "metacarpal_1_left", "metacarpal_1_right",
    "metacarpal_2_left", "metacarpal_2_right",
    "metacarpal_3_left", "metacarpal_3_right",
    "metacarpal_4_left", "metacarpal_4_right",
    "metacarpal_5_left", "metacarpal_5_right",
    "metacarpal_indet", "metacarpal_left_indet", "metacarpal_right_indet"
    

  )
),

# Level 5 -> Level 1
category_mappings_hands_5_to_1 <- list(
  hands = c(
    "capitate_left", "capitate_right",
    "carpal_indet", "carpal_left_indet", "carpal_right_indet",
    "hamate_left", "hamate_right",
    "hand_phalanx_distal_indet",
    "hand_phalanx_distal_left_1", "hand_phalanx_distal_left_2", "hand_phalanx_distal_left_3", "hand_phalanx_distal_left_4", "hand_phalanx_distal_left_5",
    "hand_phalanx_distal_left_indet",
    "hand_phalanx_distal_right_1", "hand_phalanx_distal_right_2", "hand_phalanx_distal_right_3", "hand_phalanx_distal_right_4", "hand_phalanx_distal_right_5",
    "hand_phalanx_distal_right_indet",
    "hand_phalanx_indet",
    "hand_phalanx_middle_2_left", "hand_phalanx_middle_2_right",
    "hand_phalanx_middle_3_left", "hand_phalanx_middle_3_right",
    "hand_phalanx_middle_4_left", "hand_phalanx_middle_4_right",
    "hand_phalanx_middle_5_left", "hand_phalanx_middle_5_right",
    "hand_phalanx_middle_indet",
    "hand_phalanx_middle_left_indet", "hand_phalanx_middle_right_indet",
    "hand_phalanx_proximal_1_left", "hand_phalanx_proximal_1_right",
    "hand_phalanx_proximal_2_left", "hand_phalanx_proximal_2_right",
    "hand_phalanx_proximal_3_left", "hand_phalanx_proximal_3_right",
    "hand_phalanx_proximal_4_left", "hand_phalanx_proximal_4_right",
    "hand_phalanx_proximal_5_left", "hand_phalanx_proximal_5_right",
    "hand_phalanx_proximal_indet",
    "hand_phalanx_proximal_left_indet", "hand_phalanx_proximal_right_indet",
    "hands_indet",
    "lunate_left", "lunate_right",
    "metacarpal_1_left", "metacarpal_1_right",
    "metacarpal_2_left", "metacarpal_2_right",
    "metacarpal_3_left", "metacarpal_3_right",
    "metacarpal_4_left", "metacarpal_4_right",
    "metacarpal_5_left", "metacarpal_5_right",
    "metacarpal_indet", "metacarpal_left_indet", "metacarpal_right_indet",
    "pisiform_left", "pisiform_right",
    "scaphoid_left", "scaphoid_right",
    "trapezium_left", "trapezium_right",
    "trapezoid_left", "trapezoid_right",
    "triquetrum_left", "triquetrum_right"
  )
),

# Level 3 -> Level 2
category_mappings_hands_3_to_2 <- list(
  hands_indet = c("hands_indet"),
  carpal = c("carpal_indet", "carpal_right", "carpal_left"),
  metacarpal = c("metacarpal_indet", "metacarpal_right", "metacarpal_left"),
  hand_phalanx = c("hand_phalanx_indet", "hand_phalanx_proximal", "hand_phalanx_middle", "hand_phalanx_distal")
),

# Level 4 -> Level 2
category_mappings_hands_4_to_2 <- list(
  hands_indet = c("hands_indet"),
  carpal = c(
    "carpal_indet",
    "carpal_right_indet", "scaphoid_right", "lunate_right", "triquetrum_right", "pisiform_right",
    "trapezium_right", "trapezoid_right", "capitate_right", "hamate_right",
    "carpal_left_indet", "scaphoid_left", "lunate_left", "triquetrum_left", "pisiform_left",
    "trapezium_left", "trapezoid_left", "capitate_left", "hamate_left"
  ),
  metacarpal = c(
    "metacarpal_indet",
    "metacarpal_right_indet", "metacarpal_1_right", "metacarpal_2_right", "metacarpal_3_right", "metacarpal_4_right", "metacarpal_5_right",
    "metacarpal_left_indet", "metacarpal_1_left", "metacarpal_2_left", "metacarpal_3_left", "metacarpal_4_left", "metacarpal_5_left"
  ),
  hand_phalanx = c(
    "hand_distal_left", "hand_middle_left", "hand_middle_right",
    "hand_phalanx_distal_indet", "hand_phalanx_distal_right",
    "hand_phalanx_indet", "hand_phalanx_left",
    "hand_phalanx_middle_indet", "hand_phalanx_proximal_indet",
    "hand_phalanx_right"
  )
),

# Level 5 -> Level 2
category_mappings_hands_5_to_2 <- list(
  hands_indet = c("hands_indet"),
  carpal = c(
    "carpal_indet",
    "carpal_right_indet", "scaphoid_right", "lunate_right", "triquetrum_right", "pisiform_right",
    "trapezium_right", "trapezoid_right", "capitate_right", "hamate_right",
    "carpal_left_indet", "scaphoid_left", "lunate_left", "triquetrum_left", "pisiform_left",
    "trapezium_left", "trapezoid_left", "capitate_left", "hamate_left"
  ),
  metacarpal = c(
    "metacarpal_indet",
    "metacarpal_right_indet", "metacarpal_1_right", "metacarpal_2_right", "metacarpal_3_right", "metacarpal_4_right", "metacarpal_5_right",
    "metacarpal_left_indet", "metacarpal_1_left", "metacarpal_2_left", "metacarpal_3_left", "metacarpal_4_left", "metacarpal_5_left"
  ),
  hand_phalanx = c(
    "hand_phalanx_indet",
    "hand_phalanx_proximal_indet",
    "hand_phalanx_proximal_right_indet", "hand_phalanx_proximal_1_right", "hand_phalanx_proximal_2_right", "hand_phalanx_proximal_3_right", "hand_phalanx_proximal_4_right", "hand_phalanx_proximal_5_right",
    "hand_phalanx_proximal_left_indet", "hand_phalanx_proximal_1_left", "hand_phalanx_proximal_2_left", "hand_phalanx_proximal_3_left", "hand_phalanx_proximal_4_left", "hand_phalanx_proximal_5_left",
    "hand_phalanx_middle_indet",
    "hand_phalanx_middle_right_indet", "hand_phalanx_middle_2_right", "hand_phalanx_middle_3_right", "hand_phalanx_middle_4_right", "hand_phalanx_middle_5_right",
    "hand_phalanx_middle_left_indet", "hand_phalanx_middle_2_left", "hand_phalanx_middle_3_left", "hand_phalanx_middle_4_left", "hand_phalanx_middle_5_left",
    "hand_phalanx_distal_indet",
    "hand_phalanx_distal_right_indet", "hand_phalanx_distal_right_1", "hand_phalanx_distal_right_2", "hand_phalanx_distal_right_3", "hand_phalanx_distal_right_4", "hand_phalanx_distal_right_5",
    "hand_phalanx_distal_left_indet", "hand_phalanx_distal_left_1", "hand_phalanx_distal_left_2", "hand_phalanx_distal_left_3", "hand_phalanx_distal_left_4", "hand_phalanx_distal_left_5"
  )
),

# Level 4 -> Level 3
category_mappings_hands_4_to_3 <- list(
  hands_indet = c("hands_indet"),

  carpal_indet = c("carpal_indet"),
  carpal_right = c(
    "carpal_right_indet",
    "scaphoid_right", "lunate_right", "triquetrum_right", "pisiform_right",
    "trapezium_right", "trapezoid_right", "capitate_right", "hamate_right"
  ),
  carpal_left = c(
    "carpal_left_indet",
    "scaphoid_left", "lunate_left", "triquetrum_left", "pisiform_left",
    "trapezium_left", "trapezoid_left", "capitate_left", "hamate_left"
  ),

  metacarpal_indet = c("metacarpal_indet"),
  metacarpal_right = c(
    "metacarpal_right_indet",
    "metacarpal_1_right", "metacarpal_2_right", "metacarpal_3_right", "metacarpal_4_right", "metacarpal_5_right"
  ),
  metacarpal_left = c(
    "metacarpal_left_indet",
    "metacarpal_1_left", "metacarpal_2_left", "metacarpal_3_left", "metacarpal_4_left", "metacarpal_5_left"
  ),

  hand_phalanx_indet = c("hand_phalanx_indet"),
  hand_phalanx_proximal = c("hand_phalanx_proximal_indet", "hand_phalanx_right", "hand_phalanx_left"),
  hand_phalanx_middle   = c("hand_phalanx_middle_indet", "hand_middle_right", "hand_middle_left"),
  hand_phalanx_distal   = c("hand_phalanx_distal_indet", "hand_phalanx_distal_right", "hand_distal_left")
),

# Level 5 -> Level 3
category_mappings_hands_5_to_3 <- list(
  hands_indet = c("hands_indet"),

  carpal_indet = c("carpal_indet"),
  carpal_right = c(
    "carpal_right_indet",
    "scaphoid_right", "lunate_right", "triquetrum_right", "pisiform_right",
    "trapezium_right", "trapezoid_right", "capitate_right", "hamate_right"
  ),
  carpal_left = c(
    "carpal_left_indet",
    "scaphoid_left", "lunate_left", "triquetrum_left", "pisiform_left",
    "trapezium_left", "trapezoid_left", "capitate_left", "hamate_left"
  ),

  metacarpal_indet = c("metacarpal_indet"),
  metacarpal_right = c(
    "metacarpal_right_indet",
    "metacarpal_1_right", "metacarpal_2_right", "metacarpal_3_right", "metacarpal_4_right", "metacarpal_5_right"
  ),
  metacarpal_left = c(
    "metacarpal_left_indet",
    "metacarpal_1_left", "metacarpal_2_left", "metacarpal_3_left", "metacarpal_4_left", "metacarpal_5_left"
  ),

  hand_phalanx_indet = c("hand_phalanx_indet"),
  hand_phalanx_proximal = c(
    "hand_phalanx_proximal_indet",
    "hand_phalanx_proximal_right_indet", "hand_phalanx_proximal_1_right", "hand_phalanx_proximal_2_right", "hand_phalanx_proximal_3_right", "hand_phalanx_proximal_4_right", "hand_phalanx_proximal_5_right",
    "hand_phalanx_proximal_left_indet", "hand_phalanx_proximal_1_left", "hand_phalanx_proximal_2_left", "hand_phalanx_proximal_3_left", "hand_phalanx_proximal_4_left", "hand_phalanx_proximal_5_left"
  ),
  hand_phalanx_middle = c(
    "hand_phalanx_middle_indet",
    "hand_phalanx_middle_right_indet", "hand_phalanx_middle_2_right", "hand_phalanx_middle_3_right", "hand_phalanx_middle_4_right", "hand_phalanx_middle_5_right",
    "hand_phalanx_middle_left_indet", "hand_phalanx_middle_2_left", "hand_phalanx_middle_3_left", "hand_phalanx_middle_4_left", "hand_phalanx_middle_5_left"
  ),
  hand_phalanx_distal = c(
    "hand_phalanx_distal_indet",
    "hand_phalanx_distal_right_indet", "hand_phalanx_distal_right_1", "hand_phalanx_distal_right_2", "hand_phalanx_distal_right_3", "hand_phalanx_distal_right_4", "hand_phalanx_distal_right_5",
    "hand_phalanx_distal_left_indet", "hand_phalanx_distal_left_1", "hand_phalanx_distal_left_2", "hand_phalanx_distal_left_3", "hand_phalanx_distal_left_4", "hand_phalanx_distal_left_5"
  )
),

# Level 5 -> Level 4
category_mappings_hands_5_to_4 <- list(
  hands_indet = c("hands_indet"),

  carpal_indet = c("carpal_indet"),
  carpal_right_indet = c("carpal_right_indet"),
  carpal_left_indet  = c("carpal_left_indet"),
  scaphoid_right = c("scaphoid_right"),
  lunate_right   = c("lunate_right"),
  triquetrum_right = c("triquetrum_right"),
  pisiform_right = c("pisiform_right"),
  trapezium_right = c("trapezium_right"),
  trapezoid_right = c("trapezoid_right"),
  capitate_right = c("capitate_right"),
  hamate_right   = c("hamate_right"),
  scaphoid_left = c("scaphoid_left"),
  lunate_left   = c("lunate_left"),
  triquetrum_left = c("triquetrum_left"),
  pisiform_left = c("pisiform_left"),
  trapezium_left = c("trapezium_left"),
  trapezoid_left = c("trapezoid_left"),
  capitate_left = c("capitate_left"),
  hamate_left   = c("hamate_left"),

  metacarpal_indet = c("metacarpal_indet"),
  metacarpal_right_indet = c("metacarpal_right_indet"),
  metacarpal_left_indet  = c("metacarpal_left_indet"),
  metacarpal_1_right = c("metacarpal_1_right"),
  metacarpal_2_right = c("metacarpal_2_right"),
  metacarpal_3_right = c("metacarpal_3_right"),
  metacarpal_4_right = c("metacarpal_4_right"),
  metacarpal_5_right = c("metacarpal_5_right"),
  metacarpal_1_left = c("metacarpal_1_left"),
  metacarpal_2_left = c("metacarpal_2_left"),
  metacarpal_3_left = c("metacarpal_3_left"),
  metacarpal_4_left = c("metacarpal_4_left"),
  metacarpal_5_left = c("metacarpal_5_left"),

  hand_phalanx_indet = c("hand_phalanx_indet"),
  hand_phalanx_proximal_indet = c("hand_phalanx_proximal_indet"),
  hand_phalanx_middle_indet   = c("hand_phalanx_middle_indet"),
  hand_phalanx_distal_indet   = c("hand_phalanx_distal_indet"),

  hand_phalanx_right = c(
    "hand_phalanx_proximal_right_indet",
    "hand_phalanx_proximal_1_right", "hand_phalanx_proximal_2_right", "hand_phalanx_proximal_3_right", "hand_phalanx_proximal_4_right", "hand_phalanx_proximal_5_right"
  ),
  hand_phalanx_left = c(
    "hand_phalanx_proximal_left_indet",
    "hand_phalanx_proximal_1_left", "hand_phalanx_proximal_2_left", "hand_phalanx_proximal_3_left", "hand_phalanx_proximal_4_left", "hand_phalanx_proximal_5_left"
  ),
  hand_middle_right = c(
    "hand_phalanx_middle_right_indet",
    "hand_phalanx_middle_2_right", "hand_phalanx_middle_3_right", "hand_phalanx_middle_4_right", "hand_phalanx_middle_5_right"
  ),
  hand_middle_left = c(
    "hand_phalanx_middle_left_indet",
    "hand_phalanx_middle_2_left", "hand_phalanx_middle_3_left", "hand_phalanx_middle_4_left", "hand_phalanx_middle_5_left"
  ),
  hand_phalanx_distal_right = c(
    "hand_phalanx_distal_right_indet",
    "hand_phalanx_distal_right_1", "hand_phalanx_distal_right_2", "hand_phalanx_distal_right_3", "hand_phalanx_distal_right_4", "hand_phalanx_distal_right_5"
  ),
  hand_distal_left = c(
    "hand_phalanx_distal_left_indet",
    "hand_phalanx_distal_left_1", "hand_phalanx_distal_left_2", "hand_phalanx_distal_left_3", "hand_phalanx_distal_left_4", "hand_phalanx_distal_left_5"
  )
),


#### PELVIC

# Level 2 -> Level 1
category_mappings_pelvic_2_to_1 <- list(
  pelvic = c("pelvic_indet", "sacrum", "coccyx", "os_coxae")
),

# Level 3 -> Level 1
category_mappings_pelvic_3_to_1 <- list(
  pelvic = c(
    "pelvic_indet",
    "sacrum_indet", "s1", "s2", "s3", "s4", "s5",
    "coccyx_indet", "cx1", "cx2", "cx3", "cx4",
    "os_coxae_indet", "os_coxae_right", "os_coxae_left"
  )
),

# Level 4 -> Level 1
category_mappings_pelvic_4_to_1 <- list(
  pelvic = c(
    "pelvic_indet",
    "sacrum_indet","s1" ,"s2", "s3", "s4", "s5", 
    "coccyx_indet", "cx1","cx2", "cx3", "cx4",
    "os_coxae_indet",
    "os_coxae_right_indet", "ilium_right", "ischium_right", "pubis_right",
    "os_coxae_left_indet",  "ilium_left",  "ischium_left",  "pubis_left"
  )
),

# Level 3 -> Level 2
category_mappings_pelvic_3_to_2 <- list(
  pelvic_indet = c("pelvic_indet"),
  sacrum       = c("sacrum_indet", "s1", "s2", "s3", "s4", "s5"),
  coccyx       = c("coccyx_indet", "cx1", "cx2", "cx3", "cx4"),
  os_coxae     = c("os_coxae_indet", "os_coxae_right", "os_coxae_left")
)
,
# Level 4 -> Level 2
category_mappings_pelvic_4_to_2 <- list(
  pelvic_indet = c("pelvic_indet"),
  sacrum       = c("sacrum_indet","s1" ,"s2", "s3", "s4", "s5"),
  coccyx       = c("coccyx_indet", "cx1","cx2", "cx3", "cx4"),
  os_coxae     = c(
    "os_coxae_indet",
    "os_coxae_right_indet", "ilium_right", "ischium_right", "pubis_right",
    "os_coxae_left_indet",  "ilium_left",  "ischium_left",  "pubis_left"
  )
),

# Level 4 -> Level 3
category_mappings_pelvic_4_to_3 <- list(
  pelvic_indet = c("pelvic_indet"),

  sacrum_indet = c("sacrum_indet"),
  s1 = c("s1"),
  s2 = c("s2"),
  s3 = c("s3"),
  s4 = c("s4"),
  s5 = c("s5"),

  coccyx_indet = c("coccyx_indet"),
  cx1 = c("cx1"),
  cx2 = c("cx2"),
  cx3 = c("cx3"),
  cx4 = c("cx4"),

  os_coxae_indet = c("os_coxae_indet"),
  os_coxae_right = c("os_coxae_right_indet", "ilium_right", "ischium_right", "pubis_right"),
  os_coxae_left  = c("os_coxae_left_indet",  "ilium_left",  "ischium_left",  "pubis_left")
),



#### LEGS

# Level 2 -> Level 1
category_mappings_legs_2_to_1 <- list(
  legs = c("legs_indet", "femur", "tibia", "fibula", "patella")
),

# Level 3 -> Level 1
category_mappings_legs_3_to_1 <- list(
  legs = c(
      "legs_indet", "femur_indet", "femur_right", "femur_left", "tibia_indet", "tibia_right", "tibia_left", "fibula_indet",
      "fibula_right", "fibula_left", "patella_indet", "patella_right", "patella_left"
    )
),

# Level 3 -> Level 2
category_mappings_legs_3_to_2 <- list(
  legs_indet = c("legs_indet"),
  femur = c("femur_indet", "femur_right", "femur_left"),
  tibia = c("tibia_indet", "tibia_right", "tibia_left"),
  fibula = c("fibula_indet", "fibula_right", "fibula_left"),
  patella = c("patella_indet", "patella_right", "patella_left")
),


#### FOOT

# Level 2 -> Level 1
category_mappings_foot_2_to_1 <- list(
  foot = c(
    "foot_indet", "tarsal", "metatarsal", "foot_phalanx"
  )
),

# Level 3 -> Level 1
category_mappings_foot_3_to_1 <- list(
  foot = c(
    "foot_indet", "tarsal_indet", "tarsal_right", "tarsal_left",
    "metatarsal_indet", "metatarsal_right", "metatarsal_left",
    "foot_phalanx_indet", "foot_phalanx_proximal", "foot_phalanx_middle",
    "foot_phalanx_distal"
  )
)
,
# Level 4 -> Level 1
category_mappings_foot_4_to_1 <- list(
  foot = c(
    "foot_indet", "tarsal_indet", "tarsal_right_indet", "talus_right",
    "calcaneous_right", "cuboid_right", "navicular_right",
    "lateral_cuneiform_right", "intermediate_cuneiform_right",
    "medial_cuneiform_right", "tarsal_left_indet", "talus_left",
    "calcaneous_left", "cuboid_left", "navicular_left",
    "lateral_cuneiform_left", "intermediate_cuneiform_left",
    "medial_cuneiform_left", "metatarsal_indet", "metatarsal_right_indet",
    "metatarsal_1_right", "metatarsal_2_right", "metatarsal_3_right",
    "metatarsal_4_right", "metatarsal_5_right", "metatarsal_left_indet",
    "metatarsal_1_left", "metatarsal_2_left", "metatarsal_3_left",
    "metatarsal_4_left", "metatarsal_5_left", "foot_phalanx_indet",
    "foot_phalanx_proximal_indet", "foot_proximal_right", "foot_phalanx_left",
    "foot_phalanx_middle_indet", "foot_middle_right", "foot_middle_left",
    "foot_phalanx_distal_indet", "foot_phalanx_distal_right", "foot_distal_left"
  )
),

# Level 5 -> Level 1
category_mappings_foot_5_to_1 <- list(
  foot = c(
    "foot_indet", "tarsal_indet", "tarsal_right_indet", "talus_right",
    "calcaneous_right", "cuboid_right", "navicular_right",
    "lateral_cuneiform_right", "intermediate_cuneiform_right",
    "medial_cuneiform_right", "tarsal_left_indet", "talus_left",
    "calcaneous_left", "cuboid_left", "navicular_left",
    "lateral_cuneiform_left", "intermediate_cuneiform_left",
    "medial_cuneiform_left", "metatarsal_indet", "metatarsal_right_indet",
    "metatarsal_1_right", "metatarsal_2_right", "metatarsal_3_right",
    "metatarsal_4_right", "metatarsal_5_right", "metatarsal_left_indet",
    "metatarsal_1_left", "metatarsal_2_left", "metatarsal_3_left",
    "metatarsal_4_left", "metatarsal_5_left", "foot_phalanx_indet",
    "foot_phalanx_proximal_indet", "foot_phalanx_proximal_right_indet",
    "foot_phalanx_proximal_1_right", "foot_phalanx_proximal_2_right",
    "foot_phalanx_proximal_3_right", "foot_phalanx_proximal_4_right",
    "foot_phalanx_proximal_5_right", "foot_phalanx_proximal_left_indet",
    "foot_phalanx_proximal_1_left", "foot_phalanx_proximal_2_left",
    "foot_phalanx_proximal_3_left", "foot_phalanx_proximal_4_left",
    "foot_phalanx_proximal_5_left", "foot_phalanx_middle_indet",
    "foot_phalanx_middle_right_indet", "foot_phalanx_middle_2_right",
    "foot_phalanx_middle_3_right", "foot_phalanx_middle_4_right",
    "foot_phalanx_middle_5_right", "foot_phalanx_middle_left_indet",
    "foot_phalanx_middle_2_left", "foot_phalanx_middle_3_left",
    "foot_phalanx_middle_4_left", "foot_phalanx_middle_5_left",
    "foot_phalanx_distal_indet", "foot_phalanx_distal_right_indet",
    "foot_phalanx_distal_right_1", "foot_phalanx_distal_right_2",
    "foot_phalanx_distal_right_3", "foot_phalanx_distal_right_4",
    "foot_phalanx_distal_right_5", "foot_phalanx_distal_left_indet",
    "foot_phalanx_distal_left_1", "foot_phalanx_distal_left_2",
    "foot_phalanx_distal_left_3", "foot_phalanx_distal_left_4",
    "foot_phalanx_distal_left_5"
  )
),

# Level 3 -> Level 2
category_mappings_foot_3_to_2 <- list(
  foot_indet = c(
    "foot_indet"
  ),
  tarsal = c(
    "tarsal_indet", "tarsal_right", "tarsal_left"
  ),
  metatarsal = c(
    "metatarsal_indet", "metatarsal_right", "metatarsal_left"
  ),
  foot_phalanx = c(
    "foot_phalanx_indet", "foot_phalanx_proximal", "foot_phalanx_middle",
    "foot_phalanx_distal"
  )
),

# Level 4 -> Level 2
category_mappings_foot_4_to_2 <- list(
  foot_indet = c(
    "foot_indet"
  ),
  tarsal = c(
    "tarsal_indet", "tarsal_right_indet", "talus_right", "calcaneous_right",
    "cuboid_right", "navicular_right", "lateral_cuneiform_right",
    "intermediate_cuneiform_right", "medial_cuneiform_right",
    "tarsal_left_indet", "talus_left", "calcaneous_left", "cuboid_left",
    "navicular_left", "lateral_cuneiform_left", "intermediate_cuneiform_left",
    "medial_cuneiform_left"
  ),
  metatarsal = c(
    "metatarsal_indet", "metatarsal_right_indet", "metatarsal_1_right",
    "metatarsal_2_right", "metatarsal_3_right", "metatarsal_4_right",
    "metatarsal_5_right", "metatarsal_left_indet", "metatarsal_1_left",
    "metatarsal_2_left", "metatarsal_3_left", "metatarsal_4_left",
    "metatarsal_5_left"
  ),
  foot_phalanx = c(
    "foot_phalanx_indet", "foot_phalanx_proximal_indet", "foot_proximal_right",
    "foot_phalanx_left", "foot_phalanx_middle_indet", "foot_middle_right",
    "foot_middle_left", "foot_phalanx_distal_indet", "foot_phalanx_distal_right",
    "foot_distal_left"
  )
),

# Level 5 -> Level 2
category_mappings_foot_5_to_2 <- list(
  foot_indet = c(
    "foot_indet"
  ),
  tarsal = c(
    "tarsal_indet", "tarsal_right_indet", "talus_right", "calcaneous_right",
    "cuboid_right", "navicular_right", "lateral_cuneiform_right",
    "intermediate_cuneiform_right", "medial_cuneiform_right",
    "tarsal_left_indet", "talus_left", "calcaneous_left", "cuboid_left",
    "navicular_left", "lateral_cuneiform_left", "intermediate_cuneiform_left",
    "medial_cuneiform_left"
  ),
  metatarsal = c(
    "metatarsal_indet", "metatarsal_right_indet", "metatarsal_1_right",
    "metatarsal_2_right", "metatarsal_3_right", "metatarsal_4_right",
    "metatarsal_5_right", "metatarsal_left_indet", "metatarsal_1_left",
    "metatarsal_2_left", "metatarsal_3_left", "metatarsal_4_left",
    "metatarsal_5_left"
  ),
  foot_phalanx = c(
    "foot_phalanx_indet", "foot_phalanx_proximal_indet",
    "foot_phalanx_proximal_right_indet", "foot_phalanx_proximal_1_right",
    "foot_phalanx_proximal_2_right", "foot_phalanx_proximal_3_right",
    "foot_phalanx_proximal_4_right", "foot_phalanx_proximal_5_right",
    "foot_phalanx_proximal_left_indet", "foot_phalanx_proximal_1_left",
    "foot_phalanx_proximal_2_left", "foot_phalanx_proximal_3_left",
    "foot_phalanx_proximal_4_left", "foot_phalanx_proximal_5_left",
    "foot_phalanx_middle_indet", "foot_phalanx_middle_right_indet",
    "foot_phalanx_middle_2_right", "foot_phalanx_middle_3_right",
    "foot_phalanx_middle_4_right", "foot_phalanx_middle_5_right",
    "foot_phalanx_middle_left_indet", "foot_phalanx_middle_2_left",
    "foot_phalanx_middle_3_left", "foot_phalanx_middle_4_left",
    "foot_phalanx_middle_5_left", "foot_phalanx_distal_indet",
    "foot_phalanx_distal_right_indet", "foot_phalanx_distal_right_1",
    "foot_phalanx_distal_right_2", "foot_phalanx_distal_right_3",
    "foot_phalanx_distal_right_4", "foot_phalanx_distal_right_5",
    "foot_phalanx_distal_left_indet", "foot_phalanx_distal_left_1",
    "foot_phalanx_distal_left_2", "foot_phalanx_distal_left_3",
    "foot_phalanx_distal_left_4", "foot_phalanx_distal_left_5"
  )
),

# Level 4 -> Level 3
category_mappings_foot_4_to_3 <- list(
  foot_indet = c(
    "foot_indet"
  ),
  tarsal_indet = c(
    "tarsal_indet"
  ),
  tarsal_right = c(
    "tarsal_right_indet", "talus_right", "calcaneous_right", "cuboid_right",
    "navicular_right", "lateral_cuneiform_right", "intermediate_cuneiform_right",
    "medial_cuneiform_right"
  ),
  tarsal_left = c(
    "tarsal_left_indet", "talus_left", "calcaneous_left", "cuboid_left",
    "navicular_left", "lateral_cuneiform_left", "intermediate_cuneiform_left",
    "medial_cuneiform_left"
  ),
  metatarsal_indet = c(
    "metatarsal_indet"
  ),
  metatarsal_right = c(
    "metatarsal_right_indet", "metatarsal_1_right", "metatarsal_2_right",
    "metatarsal_3_right", "metatarsal_4_right", "metatarsal_5_right"
  ),
  metatarsal_left = c(
    "metatarsal_left_indet", "metatarsal_1_left", "metatarsal_2_left",
    "metatarsal_3_left", "metatarsal_4_left", "metatarsal_5_left"
  ),
  foot_phalanx_indet = c(
    "foot_phalanx_indet"
  ),
  foot_phalanx_proximal= c(
    "foot_phalanx_proximal_indet", "foot_proximal_right", "foot_phalanx_left"
  ),
  foot_phalanx_middle = c(
    "foot_phalanx_middle_indet", "foot_middle_right", "foot_middle_left"
  ),
  foot_phalanx_distal = c(
    "foot_phalanx_distal_indet", "foot_phalanx_distal_right", "foot_distal_left"
  )
)
,
# Level 5 -> Level 3
category_mappings_foot_5_to_3 <- list(
  foot_indet = c(
    "foot_indet"
  ),
  tarsal_indet = c(
    "tarsal_indet"
  ),
  tarsal_right = c(
    "tarsal_right_indet", "talus_right", "calcaneous_right", "cuboid_right",
    "navicular_right", "lateral_cuneiform_right", "intermediate_cuneiform_right",
    "medial_cuneiform_right"
  ),
  tarsal_left = c(
    "tarsal_left_indet", "talus_left", "calcaneous_left", "cuboid_left",
    "navicular_left", "lateral_cuneiform_left", "intermediate_cuneiform_left",
    "medial_cuneiform_left"
  ),
  metatarsal_indet = c(
    "metatarsal_indet"
  ),
  metatarsal_right = c(
    "metatarsal_right_indet", "metatarsal_1_right", "metatarsal_2_right",
    "metatarsal_3_right", "metatarsal_4_right", "metatarsal_5_right"
  ),
  metatarsal_left = c(
    "metatarsal_left_indet", "metatarsal_1_left", "metatarsal_2_left",
    "metatarsal_3_left", "metatarsal_4_left", "metatarsal_5_left"
  ),
  foot_phalanx_indet = c(
    "foot_phalanx_indet"
  ),
  foot_phalanx_proximal = c(
    "foot_phalanx_proximal_indet", "foot_phalanx_proximal_right_indet",
    "foot_phalanx_proximal_1_right", "foot_phalanx_proximal_2_right",
    "foot_phalanx_proximal_3_right", "foot_phalanx_proximal_4_right",
    "foot_phalanx_proximal_5_right", "foot_phalanx_proximal_left_indet",
    "foot_phalanx_proximal_1_left", "foot_phalanx_proximal_2_left",
    "foot_phalanx_proximal_3_left", "foot_phalanx_proximal_4_left",
    "foot_phalanx_proximal_5_left"
  ),
  foot_phalanx_middle = c(
    "foot_phalanx_middle_indet", "foot_phalanx_middle_right_indet",
    "foot_phalanx_middle_2_right", "foot_phalanx_middle_3_right",
    "foot_phalanx_middle_4_right", "foot_phalanx_middle_5_right",
    "foot_phalanx_middle_left_indet", "foot_phalanx_middle_2_left",
    "foot_phalanx_middle_3_left", "foot_phalanx_middle_4_left",
    "foot_phalanx_middle_5_left"
  ),
  foot_phalanx_distal = c(
    "foot_phalanx_distal_indet", "foot_phalanx_distal_right_indet",
    "foot_phalanx_distal_right_1", "foot_phalanx_distal_right_2",
    "foot_phalanx_distal_right_3", "foot_phalanx_distal_right_4",
    "foot_phalanx_distal_right_5", "foot_phalanx_distal_left_indet",
    "foot_phalanx_distal_left_1", "foot_phalanx_distal_left_2",
    "foot_phalanx_distal_left_3", "foot_phalanx_distal_left_4",
    "foot_phalanx_distal_left_5"
  )
),

# Level 5 -> Level 4
category_mappings_foot_5_to_4 <- list(
  foot_indet = c(
    "foot_indet"
  ),
  tarsal_indet = c(
    "tarsal_indet"
  ),
  tarsal_right_indet = c(
    "tarsal_right_indet"
  ),
  talus_right = c(
    "talus_right"
  ),
  calcaneous_right = c(
    "calcaneous_right"
  ),
  cuboid_right = c(
    "cuboid_right"
  ),
  navicular_right = c(
    "navicular_right"
  ),
  lateral_cuneiform_right = c(
    "lateral_cuneiform_right"
  ),
  intermediate_cuneiform_right = c(
    "intermediate_cuneiform_right"
  ),
  medial_cuneiform_right = c(
    "medial_cuneiform_right"
  ),
  tarsal_left_indet = c(
    "tarsal_left_indet"
  ),
  talus_left = c(
    "talus_left"
  ),
  calcaneous_left = c(
    "calcaneous_left"
  ),
  cuboid_left = c(
    "cuboid_left"
  ),
  navicular_left = c(
    "navicular_left"
  ),
  lateral_cuneiform_left = c(
    "lateral_cuneiform_left"
  ),
  intermediate_cuneiform_left = c(
    "intermediate_cuneiform_left"
  ),
  medial_cuneiform_left = c(
    "medial_cuneiform_left"
  ),
  metatarsal_indet = c(
    "metatarsal_indet"
  ),
  metatarsal_right_indet = c(
    "metatarsal_right_indet"
  ),
  metatarsal_1_right = c(
    "metatarsal_1_right"
  ),
  metatarsal_2_right = c(
    "metatarsal_2_right"
  ),
  metatarsal_3_right = c(
    "metatarsal_3_right"
  ),
  metatarsal_4_right = c(
    "metatarsal_4_right"
  ),
  metatarsal_5_right = c(
    "metatarsal_5_right"
  ),
  metatarsal_left_indet = c(
    "metatarsal_left_indet"
  ),
  metatarsal_1_left = c(
    "metatarsal_1_left"
  ),
  metatarsal_2_left = c(
    "metatarsal_2_left"
  ),
  metatarsal_3_left = c(
    "metatarsal_3_left"
  ),
  metatarsal_4_left = c(
    "metatarsal_4_left"
  ),
  metatarsal_5_left = c(
    "metatarsal_5_left"
  ),
  foot_phalanx_indet = c(
    "foot_phalanx_indet"
  ),
  foot_phalanx_proximal_indet = c(
    "foot_phalanx_proximal_indet"
  ),
  foot_proximal_right = c(
    "foot_phalanx_proximal_right_indet", "foot_phalanx_proximal_1_right",
    "foot_phalanx_proximal_2_right", "foot_phalanx_proximal_3_right",
    "foot_phalanx_proximal_4_right", "foot_phalanx_proximal_5_right"
  ),
  foot_phalanx_left = c(
    "foot_phalanx_proximal_left_indet", "foot_phalanx_proximal_1_left",
    "foot_phalanx_proximal_2_left", "foot_phalanx_proximal_3_left",
    "foot_phalanx_proximal_4_left", "foot_phalanx_proximal_5_left"
  ),
  foot_phalanx_middle_indet = c(
    "foot_phalanx_middle_indet"
  ),
  foot_middle_right = c(
    "foot_phalanx_middle_right_indet", "foot_phalanx_middle_2_right",
    "foot_phalanx_middle_3_right", "foot_phalanx_middle_4_right",
    "foot_phalanx_middle_5_right"
  ),
  foot_middle_left = c(
    "foot_phalanx_middle_left_indet", "foot_phalanx_middle_2_left",
    "foot_phalanx_middle_3_left", "foot_phalanx_middle_4_left",
    "foot_phalanx_middle_5_left"
  ),
  foot_phalanx_distal_indet = c(
    "foot_phalanx_distal_indet"
  ),
  foot_phalanx_distal_right = c(
    "foot_phalanx_distal_right_indet", "foot_phalanx_distal_right_1",
    "foot_phalanx_distal_right_2", "foot_phalanx_distal_right_3",
    "foot_phalanx_distal_right_4", "foot_phalanx_distal_right_5"
  ),
  foot_distal_left = c(
    "foot_phalanx_distal_left_indet", "foot_phalanx_distal_left_1",
    "foot_phalanx_distal_left_2", "foot_phalanx_distal_left_3",
    "foot_phalanx_distal_left_4", "foot_phalanx_distal_left_5"
  )
)




  )

  ######### recal MNE ########

  category_hierarchy <- list(
    
    # Cranium 4 10 10
    "cranium_4_to_2" = category_mappings_cranium_4_to_2,
    "cranium_4_to_1" = category_mappings_cranium_4_to_1,
    "cranium_3_to_2" = category_mappings_cranium_3_to_2,
    "cranium_3_to_1" = category_mappings_cranium_3_to_1,
    "cranium_2_to_1" = category_mappings_cranium_2_to_1,
    "cranium_4_to_3" = category_mappings_cranium_4_to_3,
    "cranium_5_to_1" = category_mappings_cranium_5_to_1,
    "cranium_5_to_2" = category_mappings_cranium_5_to_2,
    "cranium_5_to_3" = category_mappings_cranium_5_to_3,
    "cranium_5_to_4" = category_mappings_cranium_5_to_4,


    # Teeth 
    "teeth_2_to_1" = category_mappings_teeth_2_to_1,
    "teeth_3_to_1" = category_mappings_teeth_3_to_1,
    "teeth_4_to_1" = category_mappings_teeth_4_to_1,
    "teeth_3_to_2" = category_mappings_teeth_3_to_2,
    "teeth_4_to_2" = category_mappings_teeth_4_to_2,
    "teeth_5_to_1" = category_mappings_teeth_5_to_1,
    "teeth_5_to_2" = category_mappings_teeth_5_to_2,
    "teeth_4_to_3" = category_mappings_teeth_4_to_3,
    "teeth_5_to_3" = category_mappings_teeth_5_to_3,
    "teeth_5_to_4" = category_mappings_teeth_5_to_4,


    #Vertebrae 3  3 3
    "vertebrae_2_to_1" = category_mappings_vertebrae_2_to_1,
    "vertebrae_3_to_1" = category_mappings_vertebrae_3_to_1,
    "vertebrae_3_to_2" = category_mappings_vertebrae_3_to_2,

    #Thorax 4  6 6
    "thorax_4_to_3" = category_mappings_thorax_4_to_3,
    "thorax_4_to_2" = category_mappings_thorax_4_to_2,
    "thorax_4_to_1" = category_mappings_thorax_4_to_1,
    "thorax_3_to_2" = category_mappings_thorax_3_to_2,
    "thorax_3_to_1" = category_mappings_thorax_3_to_1,
    "thorax_2_to_1" = category_mappings_thorax_2_to_1,

    #Shoulder 3 3 3
    "shoulder_2_to_1" = category_mappings_shoulder_2_to_1,
    "shoulder_3_to_1" = category_mappings_shoulder_3_to_1,
    "shoulder_3_to_2" = category_mappings_shoulder_3_to_2,

    #Arms 3 3 3
    "arms_2_to_1" = category_mappings_arms_2_to_1,
    "arms_3_to_1" = category_mappings_arms_3_to_1,
    "arms_3_to_2" = category_mappings_arms_3_to_2,

    #Hands 5 10 10
    "hands_2_to_1" = category_mappings_hands_2_to_1,
    "hands_3_to_1" = category_mappings_hands_3_to_1,
    "hands_4_to_1" = category_mappings_hands_4_to_1,
    "hands_3_to_2" = category_mappings_hands_3_to_2,
    "hands_4_to_2" = category_mappings_hands_4_to_2,
    "hands_5_to_1" = category_mappings_hands_5_to_1,
    "hands_5_to_2" = category_mappings_hands_5_to_2,
    "hands_4_to_3" = category_mappings_hands_4_to_3,
    "hands_5_to_3" = category_mappings_hands_5_to_3,
    "hands_5_to_4" = category_mappings_hands_5_to_4,

    #Pelvic 4 6 6
    "pelvic_2_to_1" = category_mappings_pelvic_2_to_1,
    "pelvic_3_to_1" = category_mappings_pelvic_3_to_1,
    "pelvic_3_to_2" = category_mappings_pelvic_3_to_2,
    "pelvic_4_to_1" = category_mappings_pelvic_4_to_1,
    "pelvic_4_to_2" = category_mappings_pelvic_4_to_2,
    "pelvic_4_to_3" = category_mappings_pelvic_4_to_3,

    #Leg 3 3 3
    "legs_2_to_1" = category_mappings_legs_2_to_1,
    "legs_3_to_1" = category_mappings_legs_3_to_1,
    "legs_3_to_2" = category_mappings_legs_3_to_2,

    #Foot 5 10 10
    "foot_2_to_1" = category_mappings_foot_2_to_1,
    "foot_3_to_1" = category_mappings_foot_3_to_1,
    "foot_4_to_1" = category_mappings_foot_4_to_1,
    "foot_3_to_2" = category_mappings_foot_3_to_2,
    "foot_4_to_2" = category_mappings_foot_4_to_2,
    "foot_5_to_1" = category_mappings_foot_5_to_1,
    "foot_5_to_2" = category_mappings_foot_5_to_2,
    "foot_4_to_3" = category_mappings_foot_4_to_3,
    "foot_5_to_3" = category_mappings_foot_5_to_3,
    "foot_5_to_4" = category_mappings_foot_5_to_4

  )
  

    
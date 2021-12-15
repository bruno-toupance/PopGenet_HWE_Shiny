#==============================================================================
#    ui.R : PopGenet_HWE_Shiny User-Interface
#    Copyright (C) 2021  Bruno Toupance <bruno.toupance@mnhn.fr>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#==============================================================================


library(shiny)


#==============================================================================
# shinyUI
#==============================================================================
shinyUI(

	pageWithSidebar(

#------------------------------------------------------------------------------
# Header
#------------------------------------------------------------------------------
		headerPanel("Population Genetics - HWE"), 
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Sidebar with input
#------------------------------------------------------------------------------
		sidebarPanel(

			wellPanel(
				h3("Observed counts"), 
				# numericInput(inputId = 'nb_AA', label = 'AA', value = 363),
				# numericInput(inputId = 'nb_Aa', label = 'Aa', value = 634),
				# numericInput(inputId = 'nb_aa', label = 'aa', value = 282)
				# numericInput(inputId = 'nb_AA', label = 'AA', value = 9365),
				# numericInput(inputId = 'nb_Aa', label = 'Aa', value = 2993),
				# numericInput(inputId = 'nb_aa', label = 'aa', value = 29)
				numericInput(inputId = 'nb_AA', label = 'AA', value = 96),
				numericInput(inputId = 'nb_Aa', label = 'Aa', value = 190),
				numericInput(inputId = 'nb_aa', label = 'aa', value = 143)
			),

			wellPanel(
				h3("pegas package"), 
				checkboxInput(inputId = 'pegas_perm_flag', label = 'Perform Exact Test', value = FALSE), 
				numericInput(inputId = 'pegas_B', label = 'Number of replicates', value = 1000)
			),

			wellPanel(
				h3("HardyWeinberg package"), 
				checkboxInput(inputId = 'HardyWeinberg_perm_flag', label = 'Perform Permutation Test', value = FALSE), 
				numericInput(inputId = 'HardyWeinberg_B', label = 'Number of permutations', value = 1000)
			)

		), 
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Result Panel
#------------------------------------------------------------------------------
		mainPanel(
			tabsetPanel(
				type = "tabs", 
				
#------------------------------------------------------------------------------
## Panel 'Result'
#------------------------------------------------------------------------------
				tabPanel(
					title = "Summary", 
					
					h3("Genotypes"), 
					htmlOutput("genotype_obs_table_output", container = span), 

					h3("Alleles"), 
					htmlOutput("allele_obs_table_output", container = span), 

					h3("Diversity"), 
					htmlOutput("div_table_output", container = span), 

					value = 1
				),
#------------------------------------------------------------------------------
				
#------------------------------------------------------------------------------
## Panel 'Result'
#------------------------------------------------------------------------------
				tabPanel(
					title = "Result", 
					
					h3("HWE Chisq Test"), 
					htmlOutput("hwe_table_output", container = span), 
					# htmlOutput("hwe_chisq_stat_output", container = span), 
					# htmlOutput("hwe_chisq_pvalue_output", container = span), 
					textOutput("hwe_chisq_stat_output"), 
					textOutput("hwe_chisq_pvalue_output"), 

					# hr()
					h3("Package 'pegas': HWE Exact Test"), 
					textOutput("pegas_B_output"), 
					textOutput("pegas_exact_pvalue_output"), 

					# hr()
					h3("Package 'HardyWeinberg'"), 
					htmlOutput("HardyWeinberg_output", container = span), 

					value = 1
				),
#------------------------------------------------------------------------------
				
#------------------------------------------------------------------------------
## Panel 'Plot Triangle'
#------------------------------------------------------------------------------
				tabPanel(
					title = "Plot", 
					
					splitLayout(
						plotOutput(outputId = "triangle_plot_output", height = "600px")
					), 
					
					value = 2
				)
#------------------------------------------------------------------------------
				
			)
		)
#------------------------------------------------------------------------------

	)
)


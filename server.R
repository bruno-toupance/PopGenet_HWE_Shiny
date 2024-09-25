#==============================================================================
#    server.R: PopGenet_HWE_Shiny Server
#    Copyright (C) 2024  Bruno Toupance <bruno.toupance@mnhn.fr>
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


library("shiny")

source("PopGenet_HWE.R")


#==============================================================================
# shinyServer
#==============================================================================
shinyServer(
    function(input, output, session) {
        
        # Generate population data
        pop_data <- reactive(
            {
                return(
                    compute_pop(
                        nb_AA = input$nb_AA, 
                        nb_Aa = input$nb_Aa, 
                        nb_aa = input$nb_aa
                    )
                )
            }
        )
        
        pop_pegas_data <- reactive(
            {
                return(
                    compute_pegas_pop(
                        nb_AA = input$nb_AA, 
                        nb_Aa = input$nb_Aa, 
                        nb_aa = input$nb_aa, 
                        pegas_perm_flag = input$pegas_perm_flag,
                        pegas_B = input$pegas_B
                    )
                )
            }
        )
        
        pop_HardyWeinberg_data <- reactive(
            {
                return(
                    compute_HardyWeinberg_pop(
                        nb_AA = input$nb_AA, 
                        nb_Aa = input$nb_Aa, 
                        nb_aa = input$nb_aa, 
                        HardyWeinberg_perm_flag = input$HardyWeinberg_perm_flag,
                        HardyWeinberg_B = input$HardyWeinberg_B
                    )
                )
            }
        )
        
        pop_triangle_plot_data <- reactive(
            {
                return(
                    compute_triangle_plot(
                        nb_AA = input$nb_AA, 
                        nb_Aa = input$nb_Aa, 
                        nb_aa = input$nb_aa, 
                        alpha = input$alpha
                    )
                )
            }
        )
        
        output$genotype_obs_table_output <- renderTable(
            {
                pop_data()$genotype_obs_df 
            }, 
            digits = 3
        )
        
        
        output$allele_obs_table_output <- renderTable(
            {
                pop_data()$allele_obs_df
            }, 
            digits = 3
        )
        
        
        output$div_table_output <- renderTable(
            {
                pop_data()$diversity_df
            }, 
            digits = 3
        )
        
        
        output$hwe_table_output <- renderTable(
            {
                pop_data()$hwe_df
            }, 
            digits = 3
        )
        
        output$hwe_chisq_stat_output <- renderText(
            {
                sprintf("K-squared = %.3f", pop_data()$hwe_chisq_stat)
            }
        )
        
        output$hwe_chisq_pvalue_output <- renderText(
            { 
                pvalue <- pop_data()$hwe_chisq_pvalue
                # if (pvalue < 0.001) {
                # pvalue_str <- sprintf("%.2e", pvalue)
                # } else {
                # pvalue_str <- sprintf("%.3f", pvalue)
                # }
                pvalue_str <- formatC(pvalue, digits = 3)
                sprintf("p-value = %s", pvalue_str)
            }
        )
        
        
        output$pegas_exact_pvalue_output <- renderText(
            {
                if (pop_pegas_data()$param$flag) {
                    if (pop_pegas_data()$pegas_perm_flag) {
                        pvalue <- pop_pegas_data()$pegas_data$hwe_exact_test[1, 4]
                        # if (pvalue < 0.001) {
                        # pvalue_str <- sprintf("%.2e", pvalue)
                        # } else {
                        # pvalue_str <- sprintf("%.3f", pvalue)
                        # }
                        pvalue_str <- formatC(pvalue, digits = 3)
                    } else {
                        pvalue_str <- "NA"
                    }
                    sprintf("p-value = %s", pvalue_str)
                }
            }
        )
        
        output$pegas_B_output <- renderText(
            {
                sprintf("Number of replicates = %d", pop_pegas_data()$pegas$pegas_B) 
            }
        )
        
        output$HardyWeinberg_output <- renderTable(
            {
                HW_df <- pop_HardyWeinberg_data()$HardyWeinberg_data
                HW_df[, 3] <- formatC(HW_df[, 3], digits = 3)
                HW_df
            }, 
            digits = 3
        )
        
        
        # Panel 'Plot Triangle'
        output$triangle_plot_output <- renderPlot(
            {
                new_plot <- plot_triangle(pop_triangle_plot_data())
            }
        )
        
    }
)

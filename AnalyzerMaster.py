import tkinter
import tkinter.messagebox
import customtkinter
from tkinter import filedialog
from scipy import stats
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import textwrap
import seaborn as sns

customtkinter.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        # configure window
        self.title("Analyzer Master 2")
        self.geometry(f"{1080}x{720}")
        self.tabview2 = customtkinter.CTkTabview(self, width=250)
        self.tabview2.grid(row=1, column=5, padx=(20, 0), pady=(20, 0), sticky="nsew")
        self.tabview2.add("Settings")

        # Create a frame inside the "Settings" tab
        settings_tab = self.tabview2.tab("Settings")
        settings_frame = customtkinter.CTkFrame(settings_tab)
        settings_frame.grid(row=0, column=0, padx=(20, 0), pady=(20, 0), sticky="nsew")
        entry_placeholder_text = "Control sayısı"
        entry = customtkinter.CTkEntry(settings_frame, placeholder_text=entry_placeholder_text)
        entry.grid(row=1, column=0, columnspan=2, padx=(20, 0), pady=(20, 20), sticky="nsew")

        # Create the button widget inside the settings_frame
        button_text = "Submit edin"
        button = customtkinter.CTkButton(settings_frame, fg_color="transparent", border_width=2,
                                         text_color=("gray10", "#DCE4EE"), text=button_text,
                                         command=self.start_database_search)
        button.grid(row=2, column=0, columnspan=2, padx=(20, 0), pady=(10, 10), sticky="nsew")

        # Create Positive/Negative toggle button
        self.positive_negative_var = tkinter.IntVar(value=0)  # Initial value: positive
        self.positive_negative_toggle_button = customtkinter.CTkCheckBox(master=settings_frame,
                                                                         variable=self.positive_negative_var,
                                                                         text="Negative",
                                                                         corner_radius=12)
        self.positive_negative_toggle_button.grid(row=0, column=0, padx=(40, 0), pady=(0, 0), sticky="nsew")

        # Binding the toggle button to a function
        self.positive_negative_var.trace_add("write", self.positive_negative_changed)



        # configure grid layout (4x4)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure((2, 3), weight=0)
        self.grid_rowconfigure((0, 1, 2), weight=1)

        # create sidebar frame with widgets
        self.sidebar_frame = customtkinter.CTkFrame(self, width=140, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure(4, weight=1)
        self.logo_label = customtkinter.CTkLabel(self.sidebar_frame, text="Omics analysis tools", font=customtkinter.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, column=0, padx=20, pady=(20, 10))
        self.sidebar_button_1 = customtkinter.CTkButton(self.sidebar_frame, command=self.sidebar_button_event, text="Volcano Plot")
        self.sidebar_button_1.grid(row=1, column=0, padx=20, pady=10)
        self.sidebar_button_2 = customtkinter.CTkButton(self.sidebar_frame, command=self.sidebar_button_event2, text= "Gene Ontology")
        self.sidebar_button_2.grid(row=2, column=0, padx=20, pady=10)
        self.sidebar_button_3 = customtkinter.CTkButton(self.sidebar_frame, text="Gene Set Enrichment Anaylsis", command=self.sidebar_button_event3)
        self.sidebar_button_3.grid(row=3, column=0, padx=20, pady=10)
        self.appearance_mode_label = customtkinter.CTkLabel(self.sidebar_frame, text="Appearance Mode:", anchor="w")
        self.appearance_mode_label.grid(row=5, column=0, padx=20, pady=(10, 0))
        self.appearance_mode_optionemenu = customtkinter.CTkOptionMenu(self.sidebar_frame, values=["Light", "Dark", "System"],
                                                                       command=self.change_appearance_mode_event)
        self.appearance_mode_optionemenu.grid(row=6, column=0, padx=20, pady=(10, 10))
        self.scaling_label = customtkinter.CTkLabel(self.sidebar_frame, text="UI Scaling:", anchor="w")
        self.scaling_label.grid(row=7, column=0, padx=20, pady=(10, 0))
        self.scaling_optionemenu = customtkinter.CTkOptionMenu(self.sidebar_frame, values=["80%", "90%", "100%", "110%", "120%"],
                                                               command=self.change_scaling_event)
        self.scaling_optionemenu.grid(row=8, column=0, padx=20, pady=(10, 20))

        # create main entry and button
        self.entry = customtkinter.CTkEntry(self, placeholder_text="Aramak istediğiniz geni yada GO idyi  buraya girin")
        self.entry.grid(row=3, column=1, columnspan=2, padx=(20, 0), pady=(20, 20), sticky="nsew")

        self.main_button_1 = customtkinter.CTkButton(master=self, fg_color="transparent", border_width=2,
                                                     text_color=("gray10", "#DCE4EE"),
                                                     text="Database araştırması başlat",
                                                     command=self.start_database_search)
        self.main_button_1.grid(row=3, column=3, padx=(20, 20), pady=(10, 10), sticky="nsew")

        # create textbox
        # create textbox
        self.textbox = customtkinter.CTkTextbox(self, width=450,height=300)
        self.textbox.grid(row=0, column=1, columnspan=4, padx=(20, 20), pady=(20, 0), sticky="nsew")

        # create tabview
        self.tabview = customtkinter.CTkTabview(self, width=250)
        self.tabview.grid(row=0, column=5, padx=(20, 0), pady=(20, 0), sticky="nsew")

        self.tabview.add("Mouse")
        self.tabview.add("Human")
        self.tabview.tab("Mouse").grid_columnconfigure(0, weight=1)  # configure grid of individual tabs
        self.tabview.tab("Human").grid_columnconfigure(0, weight=1)

        self.optionmenu_1 = customtkinter.CTkOptionMenu(self.tabview.tab("Mouse"), dynamic_resizing=False,
                                                        values=['GO_Biological_Process_2013', 'GO_Biological_Process_2015', 'GO_Biological_Process_2017', 'GO_Biological_Process_2017b', 'GO_Biological_Process_2018', 'GO_Biological_Process_2021', 'GO_Biological_Process_2023'])
        self.optionmenu_1.grid(row=0, column=0, padx=20, pady=(20, 10))
        self.combobox_1 = customtkinter.CTkOptionMenu(self.tabview.tab("Mouse"),
                                                    values=['GO_Cellular_Component_2013', 'GO_Cellular_Component_2015', 'GO_Cellular_Component_2017', 'GO_Cellular_Component_2017b', 'GO_Cellular_Component_2018', 'GO_Cellular_Component_2021', 'GO_Cellular_Component_2023'])

        self.combobox_2=customtkinter.CTkOptionMenu(self.tabview.tab("Mouse"),
                                                    values=['GO_Molecular_Function_2013', 'GO_Molecular_Function_2015', 'GO_Molecular_Function_2017', 'GO_Molecular_Function_2017b', 'GO_Molecular_Function_2018', 'GO_Molecular_Function_2021', 'GO_Molecular_Function_2023'])
        self.combobox_2.grid(row=2, column=0, padx=20, pady=(10, 10))
        self.combobox_1.grid(row=1, column=0, padx=20, pady=(10, 10))

        # create slider and progressbar frame
        self.slider_progressbar_frame = customtkinter.CTkFrame(self, fg_color="transparent")
        self.slider_progressbar_frame.grid(row=1, column=1, columnspan=3,padx=(20, 0), pady=(20, 0), sticky="nsew")
        self.slider_progressbar_frame.grid_columnconfigure(4, weight=1)
        self.slider_progressbar_frame.grid_rowconfigure(6, weight=1)
        self.seg_button_1 = customtkinter.CTkSegmentedButton(self.slider_progressbar_frame)
        self.seg_button_1.grid(row=0, column=0, columnspan=5 , padx=(10, 10), pady=(10, 10), sticky="ew")
        button_1 = customtkinter.CTkButton(self.seg_button_1, text="Dosya seçiniz", command=self.file_loading)
        button_1.grid(row=0, column=0, padx=(20, 10), pady=(10, 10), sticky="ew")
        self.progressbar_1 = customtkinter.CTkProgressBar(self.slider_progressbar_frame)
        self.progressbar_1.grid(row=1, column=0, columnspan=5, padx=(20, 10), pady=(10, 10), sticky="ew")
        # self.progressbar_2 = customtkinter.CTkProgressBar(self.slider_progressbar_frame)
        # self.progressbar_2.grid(row=2, column=0, padx=(20, 10), pady=(10, 10), sticky="ew")
        # self.slider_1 = customtkinter.CTkSlider(self.slider_progressbar_frame, from_=0, to=1, number_of_steps=4)
        # self.slider_1.grid(row=3, column=0, padx=(20, 10), pady=(10, 10), sticky="ew")
        self.slider_2 = customtkinter.CTkSlider(self.slider_progressbar_frame, orientation="vertical")
        self.slider_2.grid(row=0, column=5, rowspan=8, padx=(10, 100), pady=(10, 10), sticky="ns")
        self.progressbar_3 = customtkinter.CTkProgressBar(self.slider_progressbar_frame, orientation="vertical")
        self.progressbar_3.grid(row=0, column=5, rowspan=8, padx=(10, 10), pady=(10, 10), sticky="ns")

        # Add segments to the segmented button with commands


        self.appearance_mode_optionemenu.set("Dark")
        self.scaling_optionemenu.set("100%")
        self.optionmenu_1.set("Biological process")
        self.combobox_1.set("Cellular component")
        self.combobox_2.set("Molecular function")

        self.slider_2.configure(command=self.progressbar_3.set)
        self.progressbar_1.configure(mode="indeterminnate")
        self.progressbar_1.start()
        self.textbox.insert("0.0", "Welcome to Analyzer Master 2\n\n" + "\n\n")

        self.seg_button_1.set("Value 2")
        self.file_name=None
        self.volcanoplot= None
        self.geneontology= None
        self.filtered_genes= None
        self.genesetenrichment= None
        self.gene_names=None
        self.choice= None
    def positive_negative_changed(self, *args):
        # This function will be called whenever the toggle button is clicked
        # You can implement your logic here to handle positive/negative choice
        choice = self.positive_negative_var.get()
        if choice == 0:
            self.textbox.insert("end", "Up regulations selected")
            self.choice="Positive"
            # Add your positive choice logic here
        elif choice == 1:
            self.textbox.insert("end", "Down regulations selected")
            self.choice="Negative"
    def file_loading(self):
        file_name = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])  # Open file dialog
        self.file_name = file_name
        df= pd.read_csv(file_name)# Store the selected file name
        gene_names = df.iloc[:, 0]
        expression_levels = df.iloc[:, 1:7]
        expression_levels2 = df.iloc[:, 8:]

        # Perform statistical test (e.g., t-test) to calculate p-values
        t_stats, p_values = [], []
        for x, y in zip(expression_levels.values, expression_levels2.values):
            t_stat, p_value = stats.ttest_ind(x, y)
            t_stats.append(t_stat)
            p_values.append(p_value)

        # Calculate the average expression level for each gene
        average_expression = expression_levels.mean(axis=1)
        average_expression2 = expression_levels2.mean(axis=1)

        # Calculate log-fold change between two conditions
        average_expression2 = average_expression2.replace(0, np.nan).fillna(average_expression)

        log_fold_change = np.log2(average_expression2 / average_expression)

        # Create a DataFrame with the required columns
        results_df = pd.DataFrame({
            'Gene_Name': gene_names,
            'P_Value': p_values,
            'T_Statistic': t_stats,
            'Fold_Change': log_fold_change
        })

        # Filter the DataFrame to keep rows where p-value is less than 0.8
        filtered_results_df = results_df[results_df['P_Value'] > 0.8]

        # Sort the filtered DataFrame by p-value in descending order
        sorted_filtered_results_df = filtered_results_df.sort_values(by='P_Value', ascending=False)

        # Write the sorted and filtered DataFrame to a new CSV file
        self.filtered_genes= sorted_filtered_results_df
        if file_name:
            tkinter.messagebox.showinfo("Selected File", f"Selected file: {file_name}")

    def start_database_search(self):
        gene_name = self.entry.get()
        if gene_name:
            self.textbox.insert("end", f"Searching database for gene or GO id: {gene_name}\n")
            # Perform your database search using gene_name here
        else:
            tkinter.messagebox.showwarning("Input error", "Please enter a gene name.")
        self.gene_names = gene_name

    def change_appearance_mode_event(self, new_appearance_mode: str):
        customtkinter.set_appearance_mode(new_appearance_mode)

    def change_scaling_event(self, new_scaling: str):
        new_scaling_float = int(new_scaling.replace("%", "")) / 100
        customtkinter.set_widget_scaling(new_scaling_float)

    def sidebar_button_event(self):

        df = pd.read_csv(self.file_name)

        # Select the gene names and the expression levels
        gene_names = df.iloc[:, 0]
        expression_levels = df.iloc[:, 1:7]
        expression_levels2 = df.iloc[:, 8:]

        # Perform statistical test (e.g., t-test) to calculate p-values
        t_stats, p_values = [], []
        for x, y in zip(expression_levels.values, expression_levels2.values):
            t_stat, p_value = stats.ttest_ind(x, y)
            t_stats.append(t_stat)
            p_values.append(p_value)

        # Calculate the average expression level for each gene
        average_expression = expression_levels.mean(axis=1)
        average_expression2 = expression_levels2.mean(axis=1)

        # Calculate log-fold change between two conditions
        average_expression2 = np.where(average_expression2 != 0, average_expression2,average_expression)

        log_fold_change = np.log2(average_expression/ average_expression2)

        # Calculate -log10(p-value)
        minus_log_p_values = -np.log10(p_values)
        # Filter data based on log values within the range [-1.5, 1.5]
        condition = (log_fold_change >= -1.5) & (log_fold_change <= 1.5)
        filtered_log_fold_change = log_fold_change[condition]
        filtered_minus_log_p_values = minus_log_p_values[condition]
        filtered_gene_names = gene_names[condition]

        # Set significance threshold (adjust as needed)
        significance_threshold = 1.3  # -log10(0.05) = 1.3

        # Create volcano plot
        plt.figure(figsize=(12, 8))
        plt.scatter(filtered_log_fold_change, filtered_minus_log_p_values, color='blue', alpha=0.5, s=10)

        # Annotate points with gene names
        if self.gene_names == None:
            pass
        else:
            for gene, x, y in zip(filtered_gene_names, filtered_log_fold_change, filtered_minus_log_p_values):
                gene_list = self.gene_names.split(";")
                adlar = gene.upper()
                if adlar in gene_list:
                    plt.annotate(gene, (x, y), textcoords="offset points", xytext=(1, 1), ha='center')
        plt.axhline(y=significance_threshold, color='red', linestyle='--', linewidth=1.5,
                    label='Significance Threshold')
        plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--', linewidth=1, label='P-value = 0.05')
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-Log10(p-value)')
        plt.title('Volcano Plot with Gene Names (Filtered)')
        plt.legend()
        plt.grid(True)
        plt.show()

    def sidebar_button_event2(self):
        fitered_genes = self.filtered_genes
        df = fitered_genes
        # Extract the relevant columns
        genes = df.iloc[:, 0]
        fold_changes = df.iloc[:, 3]

        # Filter for negative and positive fold changes and get the gene names
        negative_genes = genes[fold_changes < 0].tolist()
        positive_genes = genes[fold_changes >= 0].tolist()

        # Perform enrichment analysis for negative fold change genes
        enr_negative = gp.enrichr(gene_list=negative_genes,  # List of genes
                                  gene_sets=['GO_Molecular_Function_2023', 'GO_Cellular_Component_2023',
                                             'GO_Biological_Process_2023'],  # List of GO term libraries
                                  organism='Mouse',  # Set organism to Mouse
                                  outdir=None)  # Do not write to disk
        # Perform enrichment analysis for positive fold change genes
        enr_positive = gp.enrichr(gene_list=positive_genes,  # List of genes
                                  gene_sets=['GO_Molecular_Function_2023', 'GO_Cellular_Component_2023',
                                             'GO_Biological_Process_2023'],  # List of GO term libraries
                                  organism='Mouse',  # Set organism to Mouse
                                  outdir=None)  # Do not write to disk

        # Initialize the dictionaries
        mydictn = {}
        mydictp = {}

        # Processing of the enrichment results for negative enrichment
        for idx, row in enr_negative.results.iterrows():
            term = row.iloc[1]
            mydictn[term] = {
                "Numberofgenes": len(row.iloc[-1]),
                "pvalues": pd.to_numeric(row.iloc[4])
            }

        # Processing the enrichment results for positive enrichment
        for idx, row in enr_positive.results.iterrows():
            term = row.iloc[1]
            mydictp[term] = {
                "Numberofgenes": len(row.iloc[-1]),
                "pvalues": pd.to_numeric(row.iloc[4])
            }

        if self.choice == "Negative":
            plt.figure(figsize=(6, 4), dpi=300)
            ax = plt.gca()
            ax.spines['top'].set_visible(True)
            ax.spines['right'].set_visible(True)

            plt.rcParams["font.size"] = 8
            plt.rcParams["font.family"] = "Times New Roman"

            down_labels = [term for term in list(mydictn.keys())[:5]]
            down_pvalues = [-np.log10(mydictn[term]["pvalues"]) for term in down_labels]
            down_ycor = [0.3 * n for n in range(1, len(down_labels) + 1)]

            plt.barh(down_ycor, down_pvalues, height=0.1, color='blue', label='Negative Enrichment')
            plt.yticks(down_ycor, [textwrap.fill(e, 30) for e in down_labels], fontsize=4)
            plt.xticks([0, 3.0, 6.0])
            plt.xlabel('-log10(p_value)')
            plt.title('Down-regulated Genes')
            plt.legend()
            plt.savefig('PRMT5_GO_plot_down.pdf', dpi=300, bbox_inches='tight')
            plt.show()

        if self.choice == "Positive":
            # Plot for up-regulated genes
            plt.figure(figsize=(3, 2), dpi=300)
            ax = plt.gca()
            ax.spines['top'].set_visible(True)
            ax.spines['right'].set_visible(True)

            plt.rcParams["font.size"] = 8
            plt.rcParams["font.family"] = "Times New Roman"

            up_labels = [term for term in list(mydictp.keys())[:5]]
            up_pvalues = [-np.log10(mydictp[term]["pvalues"]) for term in up_labels]
            up_ycor = [0.3 * n for n in range(1, len(up_labels) + 1)]

            plt.barh(up_ycor, up_pvalues, height=0.1, color='red', label='Positive Enrichment')
            plt.yticks(up_ycor, [textwrap.fill(e, 30) for e in up_labels], fontsize=4)
            plt.xticks([0, 3.0, 6.0])
            plt.xlabel('-log10(p_value)')
            plt.title('Up-regulated Genes')
            plt.legend()
            plt.savefig('PRMT5_GO_plot_up.pdf', dpi=300, bbox_inches='tight')
            plt.show()

    def sidebar_button_event3(self):
        Molecular_Function = gp.get_library(name='GO_Molecular_Function_2023', organism='Mouse')
        Biological_Process = gp.get_library(name='GO_Biological_Process_2023', organism='Mouse')
        Cellular_Component = gp.get_library(name='GO_Cellular_Component_2023', organism='Mouse')

        # Load the gene expression data
        filtred_genes = self.filtered_genes
        df = filtred_genes

        # Extract the relevant columns
        # genes = df.iloc[:, 0]
        # fold_changes = df.iloc[:, 3]
        # p_values = df.iloc[:, 1]

        # Replace zero fold changes with a small value to avoid log issues
        small_value = 1e-6
        df['Fold_Change'] = df['Fold_Change'].replace(0, small_value)

        # Calculate the ranking metric
        df['SignedLogFoldChange'] = np.sign(df['Fold_Change']) * np.log2(np.abs(df['Fold_Change']))
        df['RankMetric'] = df['SignedLogFoldChange'] * -np.log10(df['P_Value'])

        # Sort by RankMetric in descending order
        df = df.sort_values(by='RankMetric', ascending=False).reset_index(drop=True)

        # Define a function to fetch genes associated with a given GO term
        def get_genes_from_go_term(go_term):
            all_genes = {"Molecular_Function": set(), "Biological_Process": set(), "Cellular_Component": set()}
            libraries = [Molecular_Function, Biological_Process, Cellular_Component]
            for library_name, library in zip(["Molecular_Function", "Biological_Process", "Cellular_Component"],
                                             libraries):
                try:
                    # Call gp.get_library inside the function to access latest versions
                    library = gp.get_library(name=f'GO_{library_name}_2023', organism='Mouse')
                    if go_term in library:
                        all_genes[library_name] = set(library[go_term])
                except KeyError:
                    pass  # Ignore missing libraries or terms, continue searching
            return all_genes

        # Example usage
        go_term = self.gene_names

        gene_set = get_genes_from_go_term(go_term)

        # Initialize variables for computing ES
        hits = df['Gene_Name'].isin(gene_set)
        no_hits = (~hits).sum()
        es = []
        running_sum = 0

        # Compute the running sum for ES with directional contributions
        for i, hit in enumerate(hits):
            if hit:
                running_sum += df.loc[i, 'RankMetric'] / hits.sum()
            else:
                running_sum -= abs(df.loc[i, 'RankMetric']) / no_hits
            es.append(running_sum)

        ranked_list = df['RankMetric'].values
        running_enrichment_score = np.array(es)
        hit_indices = np.where(hits)[0]

        # Generate sample data (replace this with your actual data)
        ranked_genes = np.random.rand(100)
        running_enrichment_score = np.cumsum(np.random.randn(100))

        # Create the sea plot
        plt.figure(figsize=(10, 6))
        ax = sns.lineplot(x=range(1, len(ranked_genes) + 1), y=running_enrichment_score, color='blue')

        # Highlight gene hits with vertical lines
        for i, gene_hit in enumerate(ranked_genes):
            if gene_hit > 0.7:  # adjust threshold as needed
                plt.axvline(x=i + 1, color='red', linestyle='--', alpha=0.5)

        # Set labels and title
        plt.xlabel('Ranked Genes')
        plt.ylabel('Running Enrichment Score')
        plt.title('Gene Set Enrichment Analysis')

        # Add color bar
        sm = plt.cm.ScalarMappable(cmap='Blues', norm=plt.Normalize(vmin=0, vmax=1))
        sm.set_array([])  # Set an empty array since we don't have any values to map
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('Ranked Genes')

        # Show plot
        plt.grid(True)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    app = App()
    app.mainloop()
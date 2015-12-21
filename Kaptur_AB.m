%Brad Kaptur
%BE567 Final Project
%Code for Models A and B

clear

%Variables in Simulation
T_clust_mat = [];%cluster transcript
T_elem_mat = []; %active max element transcript
parsed_mat = []; %parsed mRNA transcript
AGO3_load_mat = []; %loaded AGO3
Aub_load_mat = []; %loaded Aub
Piwi_load_mat = []; %loaded Piwi
num_pp_cycles_mat = []; %number of cycles of Ping-Pong
cleaved_mat = []; %number of cleaved mRNA

%Fixed and Initial Variables
Gene = 1; %initial gene amount

%Parameters
mu_clust = 0.003; %rate of transcription cluster
mu_elem = 3; %rate of transcription active max element
delta_clust = 0.006; %rate of degradation cluster
delta_elem = 0.006; %rate of degradation active max element
k_parse = 0.003; %rate of parsing
delta_parse = 0.006; %rate of degradation of parsing
k_load_AGO3 = 0.001; %rate of loading onto AGO3
k_load_Aub = 0.001; %rate of loading onto Aub
k_load_Piwi = 0.001; %rate of loading onto Piwi
k_slicer_Aub = 0.001; %rate of slicer cleavage
k_slicer_AGO3 = 0.001; %rate of slicer cleavage
delta_cleave = 0.006; %rate of degradation of cleavage
num_parsed = 1; %when the gene is parsed, how many RNA are created

sim_length = 10000;

%many simulations
for j=1:1000
%one simulation
T_clust = 0;%cluster transcript
T_elem = 0; %active max element transcript
parsed = 0; %parsed mRNA transcript
AGO3 = 10000;
AGO3_load = 0; %loaded AGO3
Aub = 10000;
Aub_load = 0; %loaded Aub
Piwi = 10000;
cleaved = 0; %cleaved mRNA transcripts
Piwi_load = 0; %loaded Piwi
num_pp_cycles = 0; %number of cycles of Ping-Pong
time=0;


while (time < sim_length)
%one iteration of simulation    
r = rand(12,1);   

tau_T_clust = -(1/(mu_clust*Gene))*log(r(1)); %transcribe clust gene
tau_T_elem = -(1/(mu_elem*Gene))*log(r(2)); %transcribe elem gene
tau_deg_clust = -(1/(delta_clust*T_clust))*log(r(3)); %degrade clust gene
tau_deg_elem = -(1/(delta_elem*T_elem))*log(r(4)); %degrade elem gene
tau_parse = -(1/(k_parse*T_clust))*log(r(5)); %parse cluster
tau_deg_parsed = -(1/(delta_parse*parsed))*log(r(6)); %degrade parsed cluster
tau_load_Aub = -(1/(k_load_Aub*Aub*parsed))*log(r(7)); %load Aub
tau_load_Piwi = -(1/(k_load_Piwi*Piwi*parsed))*log(r(8)); %load Piwi
tau_slicer_cleave_AGO3 = -(1/(k_slicer_AGO3*T_elem*AGO3_load))*log(r(9)); %cleave at AGO3
tau_load_AGO3 = -(1/(k_load_AGO3*cleaved*AGO3))*log(r(10)); %load AGO3 with cleaved Max element
tau_slicer_cleave_Aub = -(1/(k_slicer_Aub*T_clust*Aub_load))*log(r(11)); %cleave at Aub
tau_deg_cleaved = -(1/(delta_cleave*cleaved))*log(r(12)); %degrade parsed cluster

%calculate the reaction that occurs first
tau_min = min([tau_T_clust tau_T_elem tau_deg_clust tau_deg_elem tau_parse tau_deg_parsed tau_load_Aub tau_load_Piwi tau_slicer_cleave_AGO3 tau_load_AGO3 tau_slicer_cleave_Aub tau_deg_cleaved]);

%increment time for the reaction that occurred first
time=time+tau_min;

if (time < sim_length) %only iterate if the time did not go over
    
%if tree describing list of potential reaction outcomes    
if (tau_min == tau_T_clust)
    T_clust = T_clust + 1;
elseif (tau_min == tau_T_elem)
    T_elem = T_elem + 1;
elseif (tau_min == tau_deg_clust)
    T_clust = T_clust - 1;
elseif (tau_min == tau_deg_elem)
    T_elem = T_elem - 1;
elseif (tau_min == tau_parse)
    T_clust = T_clust - 1;
    parsed = parsed + num_parsed;
elseif (tau_min == tau_deg_parsed)
    parsed = parsed - 1;
elseif (tau_min == tau_load_Aub)
    Aub_load = Aub_load + 1;
    Aub = Aub - 1;
elseif (tau_min == tau_load_Piwi)
    Piwi_load = Piwi_load + 1;
    Piwi = Piwi - 1;
elseif (tau_min == tau_slicer_cleave_AGO3)
    T_elem = T_elem - 1;
    parsed = parsed + 1;
    num_pp_cycles = num_pp_cycles + 1;
elseif (tau_min == tau_slicer_cleave_Aub)
    T_clust = T_clust - 1;
    cleaved = cleaved + 1;
elseif (tau_min == tau_load_AGO3)
    cleaved = cleaved - 1;
    AGO3_load = AGO3_load + 1;
    AGO3 = AGO3 - 1;
elseif (tau_min == tau_deg_cleaved)
    cleaved = cleaved - 1;
end

end

end

%use matrices for calculations of mean, var, etc.
T_clust_mat = [T_clust_mat T_clust];%cluster transcript
T_elem_mat = [T_elem_mat T_elem]; %active max element transcript
parsed_mat = [parsed_mat parsed]; %parsed mRNA transcript
AGO3_load_mat = [AGO3_load_mat AGO3_load]; %loaded AGO3
Aub_load_mat = [Aub_load_mat Aub_load]; %loaded Aub
Piwi_load_mat = [Piwi_load_mat Piwi_load]; %loaded Piwi
num_pp_cycles_mat = [num_pp_cycles_mat num_pp_cycles]; %number of cycles of Ping-Pong
cleaved_mat = [cleaved_mat cleaved]; %cleaved_mat

j %output counter for troubleshooting purposes
end

subplot(3,3,1)
hist(T_clust_mat)
xlabel('Number of Transcribed Cluster Genes')
ylabel('Frequency')
subplot(3,3,2)
hist(T_elem_mat)
xlabel('Number of Transcribed Element Genes')
ylabel('Frequency')
subplot(3,3,3)
hist(parsed_mat)
xlabel('Number of Parsed mRNA')
ylabel('Frequency')
subplot(3,3,4)
hist(AGO3_load_mat)
xlabel('Number of Loaded AGO3')
ylabel('Frequency')
subplot(3,3,5)
hist(Aub_load_mat)
xlabel('Number of Loaded Aub')
ylabel('Frequency')
subplot(3,3,6)
hist(Piwi_load_mat)
xlabel('Number of Loaded Piwi')
ylabel('Frequency')
subplot(3,3,7)
hist(num_pp_cycles_mat)
xlabel('Number of Ping-Pong Cycles')
ylabel('Frequency')
subplot(3,3,8)
hist(cleaved_mat)
xlabel('Number of Cleaved mRNA')
ylabel('Frequency')


%mean and variance of parameters of interest
mean_pp = mean(num_pp_cycles_mat)
var_pp = var(num_pp_cycles_mat)

mean_elem = mean(T_elem_mat)
var_elem = var(T_elem_mat)

mean_parsed = mean(parsed_mat)
var_parsed = var(parsed_mat)

mean_AGO3 = mean(AGO3_load_mat)
var_AGO3 = var(AGO3_load_mat)

mean_Aub = mean(Aub_load_mat)
var_Aub = var(Aub_load_mat)

mean_Piwi = mean(Piwi_load_mat)
var_Piwi = var(Piwi_load_mat)
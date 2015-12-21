%Brad Kaptur
%BE567 Final Project
%Code for Model D

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
AGO3_block_mat = []; %number of competitively inhibited AGO3

%Fixed and Initial Variables
Gene = 1; %initial gene amount

%Parameters
mu_clust = 0.003; %rate of transcription cluster
delta_clust = 0.006; %rate of degradation cluster
k_parse = 0.003; %rate of parsing
delta_parse = 0.006; %rate of degradation of parsing
k_load_Piwi = 0.001; %rate of loading onto Piwi
num_parsed = 1; %when the gene is parsed, how many RNA are created

sim_length = 10000;

%many simulations
for j=1:1000
%one simulation
T_clust = 0;%cluster transcript
parsed = 0; %parsed mRNA transcript
Piwi = 10000;
Piwi_load = 0; %loaded Piwi
num_pp_cycles = 0; %number of cycles of Ping-Pong
time=0;


while (time < sim_length)
%one iteration of simulation    
r = rand(13,1);   

tau_T_clust = -(1/(mu_clust*Gene))*log(r(1)); %transcribe clust gene
tau_deg_clust = -(1/(delta_clust*T_clust))*log(r(3)); %degrade clust gene
tau_parse = -(1/(k_parse*T_clust))*log(r(5)); %parse cluster
tau_deg_parsed = -(1/(delta_parse*parsed))*log(r(6)); %degrade parsed cluster
tau_load_Piwi = -(1/(k_load_Piwi*Piwi*parsed))*log(r(8)); %load Piwi

tau_min = min([tau_T_clust tau_deg_clust tau_parse tau_deg_parsed tau_load_Piwi]);

time=time+tau_min;

if (time < sim_length) %only iterate if the time did not go over
    
if (tau_min == tau_T_clust)
    T_clust = T_clust + 1;
elseif (tau_min == tau_deg_clust)
    T_clust = T_clust - 1;
elseif (tau_min == tau_parse)
    T_clust = T_clust - 1;
    parsed = parsed + num_parsed;
elseif (tau_min == tau_deg_parsed)
    parsed = parsed - 1;
elseif (tau_min == tau_load_Piwi)
    Piwi_load = Piwi_load + 1;
    Piwi = Piwi - 1;
end

end

end
T_clust_mat = [T_clust_mat T_clust];%cluster transcript
parsed_mat = [parsed_mat parsed]; %parsed mRNA transcript
Piwi_load_mat = [Piwi_load_mat Piwi_load]; %loaded Piwi
num_pp_cycles_mat = [num_pp_cycles_mat num_pp_cycles]; %number of cycles of Ping-Pong

j
end

subplot(1,3,1)
hist(T_clust_mat)
xlabel('Number of Transcribed Cluster Genes')
ylabel('Frequency')
subplot(1,3,2)
hist(parsed_mat)
xlabel('Number of Parsed mRNA')
ylabel('Frequency')
subplot(1,3,3)
hist(Piwi_load_mat)
xlabel('Number of Loaded Piwi')
ylabel('Frequency')



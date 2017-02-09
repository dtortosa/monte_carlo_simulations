monte_carlo = function(ini_freq_A, n.indv, n.simul, n.gens, sA, sB, output_name){
    #Monte carlo simulation to model long-term allele frequency changes in a autosomal locus with two alleles in a finite population of diploid organisms with sexual reproduction.  

    ini_freq_B = 1-ini_freq_A #calculate the initial frequency of allele B as the differences between 1 and the initial frequency of allele A. 
    
    all_freq_A = matrix(0, nrow=length(n.gens)+1, ncol=length(n.simul))
    all_freq_B = matrix(0, nrow=length(n.gens)+1, ncol=length(n.simul)) 

    for (i in n.simul){ #loop for each simulation (run as many times as simulations we have established in n.simul parameter). 

    copies_A = round((n.indv*2)*ini_freq_A)
    copies_B = round((n.indv*2)*ini_freq_B) #calculate the copy numbers of each allele by means of multiply its frequency by the total of alleles (double of total individuals)

    #create a population data set with cbind
    pop = cbind(seq(1:n.indv), #id of each individual
        rep(x=NA, times=n.indv), 
        rep(x=NA, times=n.indv)) #one column for cromosome. These column will have the alleles. 
    pop = as.data.frame(pop) 
    names(pop) <- c("id", "cromosome_1", "cromosome_2") #convert to data.frame and give name to columns. 

    #calculate the number of copies of each allele in cromosome 1 according to allele frequency parameters
    if (copies_A>nrow(pop)){ #If number of copies of A allele is higher than the number of indivduals:
        n_copA_c1=nrow(pop) #All individual will have A allele in at least one cromosome. 
    } else { #If not:
        if(copies_B>nrow(pop)){ #If number of copies of B allele is higher than the number of indivduals: 
            n_copA_c1=0 #No one individual will have A allele in one of the cromosomes           
        } else { #If not, it is to say, A and B have the same number of copies:
            n_copA_c1 = sample(x=1:copies_A, size=1) #number of copies of A will be a random number selected from a vector from 1 to the number of A copies. It could be selected 1, 2..... number of copies of A.                 
        }

    }

    #create columns with allele copies for both cromosomes
    #"0" will be A allele, whilst "1" will be B allele.
    #cromosome_1
    copies_A_cromosome_1 = rep(x=0, times=n_copA_c1) #create a vector with many zeros as copies of A we have in the cromosome 1 (zero=A). 
    copies_B_cromosome_1 = rep(x=1, times=nrow(pop)-n_copA_c1) #The same with allele B, the number of B copies is calculated by means of the difference between number of individuals and the number of copies of A in cromosome 1. 
    cromosome_index_first_copies = sample(2:3, 1) #cromosome that will have these copies. It will be selected random from the two cromosomes (columns 2 and 3, because the row 1 is "id")
    pop[,cromosome_index_first_copies] = sample(c(copies_A_cromosome_1, copies_B_cromosome_1)) #join and mix these vectors. And are included in the column "cromosome_index_first_copies" of the data frame.
    #cromosome 2
    copies_A_cromosome_2 = rep(x=0, times=copies_A - length(copies_A_cromosome_1)) #the number of copies of A in cromosome 2 is the difference between the total number of copies of A and the copies in the cromosome 1. 
    copies_B_cromosome_2 = rep(x=1, times=nrow(pop)-length(copies_A_cromosome_2)) #Like in the cromosome 1, we calculated the number of copies of B by means of the difference of total individuals and the copies of A in the cromosome 2.
    cromosome_index_second_copies = 2:3
    cromosome_index_second_copies = cromosome_index_second_copies[ !cromosome_index_second_copies==cromosome_index_first_copies] #select the restant cromosome, which does not have copies.   
    pop[,cromosome_index_second_copies] = sample(c(copies_A_cromosome_2, copies_B_cromosome_2)) #Column of cromosome "cromosome_index_second_copies" will be both vectors joined and mixed. 

    #Calculate percentages of selection for each genotype
    wAA = 1-sA
    wAB = 1-sA
    wBB = 1-sB

    percent_indv_drop_AA = (1-wAA)*100
    percent_indv_drop_AB = (1-wAB)*100
    percent_indv_drop_BB = (1-wBB)*100

    #Finally, before to begin the simulation we are going to establish the initial frequencies of A and B
    freq_A = ini_freq_A
    freq_B = ini_freq_B

    for (j in n.gens){ #loop for each generation (iside each simulation): 
        mat_vector = 1:n.indv #the vector for select the random pairs, length similar to N of population
        c1 = NULL 
        c2 = NULL #create empty vectors for each cromosome 

        for (k in 1:n.indv){ #as many times as individuals exist in the population
            pair = sample(1:n.indv, 2) #create new matings,as many as the number of deaths. Random selection of the individuals that will be reproductives
            c1 = append(c1, pop[pair,][1,sample(2:3,1)]) #Select the first partner: Select two rows from pop usind pair as index. Of these rows, select the first and take the value of the two or third column by means of the result of sample (like flip a coin)
            c2 = append(c2, pop[pair,][2,sample(2:3,1)]) #Select the second partnet: the same for the second partner but taking the second row form the selection of pop usind pair. 
            
        }

        new.gen = cbind(rep(1:n.indv), c1,c2) #join A and B, along with a vector with a length similar to the number of individuals. 
        new.gen = as.data.frame(new.gen) #conert in a data.frame
        names(new.gen) <- c("id", "cromosome_1", "cromosome_2") #give the numbers ot the columns

        #Drop individuals dead by selection
        new.gen_AA = new.gen[new.gen$cromosome_1=="0" & new.gen$cromosome_2=="0",]
        new.gen_AB = new.gen[new.gen$cromosome_1=="0" & new.gen$cromosome_2=="1" | new.gen$cromosome_1=="1" & new.gen$cromosome_2=="0" ,]
        new.gen_BB = new.gen[new.gen$cromosome_1=="1" & new.gen$cromosome_2=="1",]

        if (nrow(new.gen_AA)>0 & percent_indv_drop_AA>0){
            id_indv_drop_AA = new.gen_AA[sample(1:nrow(new.gen_AA), round((percent_indv_drop_AA*nrow(new.gen_AA))/100)),]$id
        } else {
            id_indv_drop_AA = 0
        }

        if (nrow(new.gen_AB)>0 & percent_indv_drop_AB>0){
            id_indv_drop_AB = new.gen_AB[sample(1:nrow(new.gen_AB), round((percent_indv_drop_AB*nrow(new.gen_AB))/100)),]$id
        } else {
            id_indv_drop_AB = 0
        }

        if (nrow(new.gen_BB)>0 & percent_indv_drop_BB>0){
            id_indv_drop_BB = new.gen_BB[sample(1:nrow(new.gen_BB), round((percent_indv_drop_BB*nrow(new.gen_BB))/100)),]$id
        } else {
            id_indv_drop_BB = 0
        }

        n_indv_drop_total = nrow(new.gen[(new.gen$id %in% id_indv_drop_AA | new.gen$id %in% id_indv_drop_AB | new.gen$id %in% id_indv_drop_BB),])

        new.gen = new.gen[!(new.gen$id %in% id_indv_drop_AA | new.gen$id %in% id_indv_drop_AB | new.gen$id %in% id_indv_drop_BB),]

        #Create new individuals from the parental generation in the same number of deads, by means of creation of new random matings. It doesn't matter that pairs are repetead, we want as many matings as deaths. The important point is that the matings have to be random. 
        c1_news = NULL 
        c2_news = NULL #create each cromosome of the news births  

        for (k in 1:n_indv_drop_total){ #for each individual dropped we created other            
            pair = sample(1:n.indv, 2) #create new matings,as many as the number of deaths. Random selection of the individuals that will be reproductives
            c1_news = append(c1_news, pop[pair,][1,sample(2:3,1)]) #Select the first partner: Select two rows from pop usind pair as index. Of these rows, select the first and take the value of the two or third column by means of the result of sample (like flip a coin)
            c2_news = append(c2_news, pop[pair,][2,sample(2:3,1)]) #Select the second partnet: the same for the second partner but taking the second row form the selection of pop usind pair. 
        }

        #create a vector with id of all indivudals killed
        id_selection = c(id_indv_drop_AA, id_indv_drop_AB, id_indv_drop_BB)
        
        if (length(id_selection[!id_selection==0])>0){#if the number of death is higher than 0:
            id_selection = id_selection[!id_selection==0]
            news = cbind(id_selection, c1_news, c2_news)
            colnames(news) = c("id", "cromosome_1", "cromosome_2")
            new.gen = rbind(new.gen, news)
            new.gen = new.gen[sample(nrow(new.gen)),]
        } else { #if not we don't have to change anything
            new.gen = new.gen
        }

        #Calculate frequencies
        freq_A = append(freq_A, (nrow(new.gen[new.gen$cromosome_1==0,]) + nrow(new.gen[new.gen$cromosome_2==0,]))/(nrow(new.gen)*2)) #calculate the frequency of allele A as the number of rows with 0 dividided by (the total of row multiplied by 2; totals alleles)
        freq_B = append(freq_B, (nrow(new.gen[new.gen$cromosome_1==1,]) + nrow(new.gen[new.gen$cromosome_2==1,]))/(nrow(new.gen)*2)) #The same for B but using its index (1)
        pop = new.gen #the new generation replace the last generation. 
    }

    #sabe the frequencies of all generations of 1 simulation in two data.frames, one for each allele. As many rows as generations, and as many columns as simulations
    all_freq_A[,i] = freq_A
    all_freq_B[,i] = freq_B

    }

    #write the the frequencies of each allele
    write.csv(all_freq_A, paste("freq_A", output_name, sep="_"), row.names=FALSE)
    write.csv(all_freq_B, paste("freq_B", output_name, sep="_"), row.names=FALSE)

    #Calculate time to fixation. 
    time_to_fix_A = NULL
    for (i in 1:ncol(all_freq_A)){ #for each simulation (1 simulation for column)
        if (!1 %in% all_freq_A[,i]){ #if there is not a frequency value of 1
            time_to_fix_A =  append(time_to_fix_A, NA)
        } else { #if not, it is to say, there is fijxation, include the rows with the first 1 (generation of fixation)
            time_to_fix_A = append(time_to_fix_A, min(which(all_freq_A[,i]==1)))
        }
    }

    time_to_fix_B = NULL #The same with allele B
    for (i in 1:ncol(all_freq_B)){
        if (!1 %in% all_freq_B[,i]){
            time_to_fix_B =  append(time_to_fix_B, NA)
        } else { 
            time_to_fix_B = append(time_to_fix_B, min(which(all_freq_B[,i]==1)))
        }
    }

    time_to_fix = as.data.frame(cbind(time_to_fix_A,time_to_fix_B))
    names(time_to_fix)[1] = "allele_A"
    names(time_to_fix)[2] = "allele_B"

    #write in a csv
    write.csv(time_to_fix, paste("time_to_fix", output_name, sep="_"), row.names=FALSE)

    #Calucalte mean and sd for time to fixation
    time_to_fix_mean_A = sum(na.omit(time_to_fix_A))/length(na.omit(time_to_fix_A)) 
    time_to_fix_sd_A = sd(na.omit(time_to_fix_A))
    time_to_fix_mean_B = sum(na.omit(time_to_fix_B))/length(na.omit(time_to_fix_B)) 
    time_to_fix_sd_B = sd(na.omit(time_to_fix_B))

    #Return in the terminal a message with time to fixation, and number of generation in which fixation occur for each allele. 
    if (length(na.omit(time_to_fix_A) > 0)){
        print(paste(paste("A was fixed in", paste(length(na.omit(time_to_fix_A)), "simulations." ,sep=" "), sep=" "), paste("The number of generations to fixation was", paste(round(time_to_fix_mean_A, digits=2), paste("\u00b1", round(time_to_fix_sd_A, digits = 2), sep=" "), sep=" "), sep=" "), sep=" ")) 
    } else {
        print(paste("A was fixed in", paste(length(na.omit(time_to_fix_A)), "simulations." ,sep=" "), sep=" ")) 
    }

    if (length(na.omit(time_to_fix_B) > 0)){
        print(paste(paste("B was fixed in", paste(length(na.omit(time_to_fix_B)), "simulations." ,sep=" "), sep=" "), paste("The number of generations to fixation was", paste(round(time_to_fix_mean_B, digits=2), paste("\u00b1", round(time_to_fix_sd_B, digits = 2), sep=" "), sep=" "), sep=" "), sep=" ")) 
    } else {
        print(paste("B was fixed in", paste(length(na.omit(time_to_fix_B)), "simulations." ,sep=" "), sep=" ")) 
    }
    
}





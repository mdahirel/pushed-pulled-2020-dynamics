########TO CLEAN UP 


###solution retenue: 

###in dication que le trt reduit connectivite: correl G0 fec disp
### % sent way increases with neggs in ps, but not pl, indicate no impediment no mvt in pl

### indication that there is DDD anyway, independent of trt, and that both are pushed, just differentlu
### proba of advanced as f is large (use of is large = is > 100 because N at front bimodal, cf suppl)


test1<-data_front %>% group_by(landscapeID) %>% arrange(Generation) %>% mutate(front_next=lead(front))

testx <- data_popsize %>% 
  select(landscapeID,Macro,Generation,Neggs_estimated,Patch)

testy=left_join(test1,testx) %>% 
  filter(Patch==front & Macro !="MacroG13") %>% 
  mutate(advance=1+1*(front_next>=front)+1*(front_next>front)+1*(front_next>(front+1))) %>% 
  mutate(NeggsSD=Neggs_estimated/sd(Neggs_estimated,na.rm=TRUE)) %>% 
  mutate(is.large=Neggs_estimated>100) %>% 
  mutate(N_cat=round(Neggs_estimated/22.5)+1) %>% 
  ungroup() %>% 
  mutate(location=front_next-front) %>% mutate(location=location-min(location,na.rm=TRUE))


mod=brm(advance>2~Dynamics*is.large+(is.large|landscapeID),
        data=filter(testy,Macro=="MacroG1"),family=bernoulli,chains=2,iter=200,prior=c(set_prior("normal(0,1.5)",class=c("sd","b"))))

mod=brm(advance>2~Dynamics*mo(N_cat)+(mo(N_cat)|landscapeID),
        data=filter(testy,Macro=="MacroG1"),family=bernoulli,chains=2,iter=200,prior=c(set_prior("normal(0,1.5)",class=c("sd","b"))))



testz=testy %>% group_by(landscapeID,Generation,Dynamics) %>% 
  summarise(advance=mean(advance,na.rm=TRUE),Patch=mean(Patch,na.rm=TRUE),Neggs=mean((Neggs_estimated),na.rm=TRUE),
            NeggsSE=plotrix::std.error(log(Neggs_estimated),na.rm=TRUE)) %>% 
  mutate(NeggsSE=NeggsSE+0.01)

testy %>% filter(Generation<13) %>% 
  ggplot()+geom_boxplot(aes(x=factor(round(Neggs_estimated/450,1)),y=as.numeric(advance>2),col=Dynamics))+facet_wrap(~Dynamics+Macro)


### be careful because failed advances are repeated in the data set (same patch is the front several successive gens)
### also look at whether it's only the 1 front patch we should look or the 2-3 (based on front width from profile model)

mod=brm(advance>2~Dynamics*logit(Neggs_estimated/450)+(logit(Neggs_estimated/450)|landscapeID)+(1|landscapeID:Patch),
        data=testy,family=bernoulli,chains=2,iter=200,prior=c(set_prior("normal(0,1.5)",class=c("sd","b"))))



mod=brm(advance>2~Dynamics+s(Neggs,by=Dynamics)+(1|landscapeID),
        data=testz,family=bernoulli,chains=2,iter=200,prior=c(set_prior("normal(0,1.5)",class=c("sd","b"))))

expand_grid(Neggs=c(1:450),Dynamics=unique(testz$Dynamics),landscapeID=testz$landscapeID[1]) %>% 
  add_fitted_draws(mod=mod,re_formula = NA) %>% 
  compare_levels(variable=.value,by=Dynamics) %>% 
  ggplot() + 
  stat_lineribbon(aes(x=Neggs,y=.value))


mod=brm(bf(Neggs_estimated|trials(450)~1,theta1~Dynamics),
        data=testy,family=mixture(binomial,binomial),chains=2,iter=200,prior=c(set_prior("normal(0,1)",dpar=c("theta1"))))


##TO DO x=pop time t patch i, y = time for patch i+1 to reach >100 individuals

##
##
##

testa <- data_popsize %>% 
  select(landscapeID,Macro,Generation,N_t=Neggs_estimated,Patch)

testaa <- data_popsize %>% 
  select(landscapeID,Macro,Generation,N_t1=Neggs_estimated,Patch) %>% 
  mutate(Generation=Generation-1,Patch=Patch-1)


testb <- data_popsize %>% 
  select(landscapeID,Macro,Generation,Patch,Neggs_estimated) %>% 
  filter(Neggs_estimated>=22.5) %>%   ##5% of hosts 
  mutate(Patch=Patch-1) %>% 
  group_by(landscapeID,Macro,Patch) %>% 
  summarise(timeto100=min(Generation))


tab=left_join(test1,testa) %>% 
  left_join(testb) %>% 
  filter(Patch==front & Macro !="MacroG13" & Generation<13) %>% 
  mutate(is.large = factor(N_t >100)) %>% 
  mutate(timeto100=replace_na(timeto100,14)) %>% 
  mutate(dtimeto100=timeto100-Generation) %>% 
  filter(dtimeto100>0)   ##front qui ont reculé, trouver un moyen de les réintegrer plus tard

  
mod=brm(dtimeto100-1|cens(timeto100==14)~Dynamics*is.large+(1|landscapeID),
        data=filter(tab,Macro=="MacroG1"),family=poisson,chains=2,iter=200,prior=c(set_prior("normal(0,1)",class=c("sd","b","Intercept"))))

mod=brm(dtimeto100==1~is.pushed*logit(N_t/450)+(logit(N_t/450)|landscapeID)+(1|landscapeID:Generation),
        data=tab,family=bernoulli,chains=2,iter=200,prior=c(set_prior("normal(0,1)",class=c("sd","b","Intercept"))))


###delay to colonisation: 0 = no delay new patch colonised next, gen; 1  1 gen delay...



testa <- data_popsize %>% 
  select(landscapeID,Macro,Generation,N_t=Neggs_estimated,Patch)

testaa <- data_popsize %>% 
  select(landscapeID,Macro,Generation,N_t1=Neggs_estimated,Patch) %>% 
  mutate(Generation=Generation-1,Patch=Patch-1)

tab=left_join(test1,testa) %>% 
  left_join(testaa) %>% 
  filter(Patch==front & Macro !="MacroG13" & Generation<13) %>% 
  mutate(is.large = factor(N_t >100))

mod=brm(bf(N_t1|trials(450)~(1|landscapeID),theta1~Dynamics*is.large+(1|landscapeID),zi1~Dynamics*is.large+(1|landscapeID)),
        data=filter(tab,Macro=="MacroG1"),
        family=mixture(zero_inflated_binomial,binomial),
        chains=2,iter=500,prior=c(set_prior("normal(0,1.5)",dpar=c("theta1","zi1"))))



mod=brm(bf(N_t1|trials(450)~1,
           family=mixture(zero_inflated_binomial,zero_inflated_binomial),
           theta1~Dynamics*is.large,
           nlf(zi1~a),nlf(zi2~a),nlf(a~inv_logit(logita)),
           logita~Dynamics*is.large,nl=TRUE
           ),
        data=filter(tab,Macro=="MacroG1"),
        chains=2,iter=200,
        prior=c(
          set_prior("normal(0,1.5)",dpar="theta1"),
          set_prior("normal(0,1.5)",nlpar="logita")
          )
        )



mod=brm(bf(N_t1|trials(450)~Dynamics*logit(N_t/450),family=binomial),
        data=filter(tab,Macro=="MacroG1" & N_t1>0),
        chains=2,iter=500,
        prior=c(
          set_prior("normal(0,1.5)",class="b")
        )
)



mod=brm(bf(N_t1/N_t~Dynamics,family=lognormal),
        data=filter(tab,Macro=="MacroG1" & N_t1>0),
        chains=2,iter=500,
        prior=c(
          set_prior("normal(0,1.5)",class="b")
        )
)




data_G0 <- filter(data_popsize,Generation==1 & Patch<=4) %>% 
  select(landscapeID,Mix,Dynamics,is.pushed,Patch,Neggs_estimated,Macro) %>% 
  ungroup() %>% 
  pivot_wider(.,names_from="Patch",values_from=c("Neggs_estimated"),names_prefix = "P") %>% 
  filter(is.na(P0)==FALSE) %>%   ##the one macro with missing info for residents
  mutate(P1=replace_na(P1,0)) %>% 
  mutate(P2=replace_na(P2,0)) %>% 
  mutate(P3=replace_na(P3,0)) %>% 
  mutate(P4=replace_na(P4,0))  

data_G0$y <- with(data_G0, cbind(P0,P1,P2,P3,P4))
data_G0$tot=data_G0 %>% select(P0,P1,P2,P3,P4) %>% rowSums()

mod=brm(y|trials(tot)~Dynamics,data=data_G0,family=multinomial,iter=200,chains=2)



mod=brm(y|trials(tot)~is.pushed+(is.pushed|Macro),data=data_G0,family=multinomial,iter=2000,chains=4,
        prior=c(
          set_prior("normal(0,1)",class="sd",dpar=c("muP1","muP2","muP3","muP4")),
          set_prior("normal(0,1.5)",class="b",dpar=c("muP1","muP2","muP3","muP4"))
        )
)

plot(conditional_effects(mod,categorical=TRUE))


data_G0 <- filter(data_popsize,Generation==1 & Patch<=4) %>% 
  select(landscapeID,Mix,Dynamics,is.pushed,Patch,Neggs_estimated,Macro) %>% 
  ungroup() %>% 
  pivot_wider(.,names_from="Patch",values_from=c("Neggs_estimated"),names_prefix = "P") %>% 
  filter(is.na(P0)==FALSE) %>%   ##the one macro with missing info for residents
  pivot_longer(cols=c(P0,P1,P2,P3,P4)) %>% 
  mutate(value=replace_na(value,0)) %>% 
  group_by(landscapeID,Macro) %>% 
  arrange(name) %>% 
  mutate(cumN=cumsum(value),totN=sum(value)) %>% 
  mutate(Patch=as.numeric(factor(name)))
  
test=data_G0 %>% uncount(weights=value)

mod=brm(cumN|trials(totN)~mo(Patch)*Dynamics,data=data_G0,family=binomial,iter=100,chains=2)


mod=brm(cumN|trials(totN)~mo(Patch)*Dynamics+(1|landscapeID),data=data_G0,family=binomial,iter=100,chains=2)


mod=brm(Patch~Dynamics+(1|landscapeID),data=test,family=cumulative,iter=100,chains=2)

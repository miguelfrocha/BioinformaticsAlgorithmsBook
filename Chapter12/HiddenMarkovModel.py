# Some of the methods on this class adapt code from:
# https://github.com/jason2506/PythonHMM (Retrieved October 2017)
# https://github.com/ananthpn/pyhmm (Retrieved October 2017)

class HiddenMarkovModel:
    
    def __init__(self, init_probs, emission_probs, trans_probs):
        """Create a constructor based on three different attributes: probability of start states; emission probabilities matrix; transition probabilities matrix.
        Both emission and transition probability matrices can be implemented as dictionaries of dictionaries. States and symbols are represented as lists and can be infered from the probabilities.
        """
        self.initstate_probs = init_probs
        self.emission_probs = emission_probs
        self.transition_probs = trans_probs
        self.states = self.emission_probs.keys()
        self.symbols = list(self.emission_probs[list(self.emission_probs.keys())[0]].keys())

    
    def get_init_prob(self, state):
        '''Get initial probability of a given state'''
        if state in self.states:
            return (self.initstate_probs[state])
        else:
            return 0

    def get_emission_prob(self, state, symbol):
        '''Get probability of a given state to emit a symbol'''
        if state in self.states and symbol in self.symbols:
            return (self.emission_probs[state][symbol])
        else:
            return 0
        
    def get_transition_prob(self, state_orig, state_dest):
        '''Get probability of transition from a origin state to destination state'''
        if state_orig in self.states and state_dest in self.states:
            return (self.transition_probs[state_orig][state_dest])
        else:
            return 0
    
    
    def set_init_prob(self, state, p):
        '''Set initial probability of a given state'''
        if state in self.states:
            self.initstate_probs[state] = p

    def set_emission_prob(self, state, symbol, p):
        '''Set probability of a given state emit a symbol'''
        if state in self.states and symbol in self.symbols:
            self.emission_probs[state][symbol] = p

    def set_transition_prob(self, state_orig, state_dest, p):
        '''Set probability of transition from a origin state to destination state'''
        if state_orig in self.states and state_dest in self.states:
            self.transition_probs[state_orig][state_dest] = p
        
    def joint_probability(self, sequence, path):
        '''Given an observed sequence and a corresponding state path calculate the probability of the sequence given the path under the model'''
        seq_len = len(sequence)
        if seq_len == 0:
            return None
        
        path_len = len(path)
        if seq_len != path_len:
            print ("Observed sequence and state path of different lenghts!")
            return None
        
        prob = self.get_init_prob(path[0]) * self.get_emission_prob(path[0], sequence[0])
        for i in range(1, len(sequence)):
            prob = prob * self.get_transition_prob(path[i-1], path[i]) * self.get_emission_prob(path[i], sequence[i])
        
        return prob
    
    
    
    def forward(self, sequence):
        '''Given an observed sequence calculate the list of forward probabilities of the sequence
        using the chain rules'''
        seq_len = len(sequence)
        if seq_len == 0:
            return []

        # calculate the product of the initial probability of each state and the first symbol of the sequence
        prob_list = [{}]
        for state in self.states:
            prob_list[0][state] = self.get_init_prob(state) * self.get_emission_prob(state, sequence[0])
        # iterate through the sequence and for each state multiply by the transition probability with any other of the possibles states; this corresponds to a jump to a new state
        # once in this new state multiply by the corresponding emission probability of the sequence symbol in that state
        for i in range(1, seq_len):
            prob_list.append({})
            for state_dest in self.states:
                prob = 0
                for state_orig in self.states:
                    prob += prob_list[i-1][state_orig] * self.get_transition_prob(state_orig, state_dest)
                prob_list[i][state_dest] = prob * self.get_emission_prob(state_dest, sequence[i])
                
        return prob_list
        # in alternative one can return the sum of all probabilities for the last symbol of the sequence
        #return sum(probs_list[-1].values())
    
    
    def backward(self, sequence):
        '''Given an observed sequence calculate the list of backward probabilities of the sequence
        This is the probability of starting in a state si at position t of the sequence
        and generate the remainder of the sequence from t to end. For this we start from the end of the sequence and compute as a variant of the forward algorithm.
        The recurrence formula for a position i and state k is given by the sum of the product of transition probability of any state to stake k times the emission probability the observed sequence symbol at position i being emitted by any of the states times backward probability of the state at position i+1. Probability of the states at the last sequence symbol is 1.
        '''
        seq_len = len(sequence)
        if seq_len == 0:
            return []

        beta = [{}]
        for state in self.states:
            beta[0][state] = 1

        for i in range(seq_len - 1, 0, -1):
            beta.insert(0, {})
            for state_orig in self.states:
                prob = 0
                for state_dest in self.states:
                    prob += beta[1][state_dest] * self.get_transition_prob(state_orig, state_dest) * self.get_emission_prob(state_dest, sequence[i])  
                beta[0][state_orig] = prob
        return beta
    
    
    def viterbi(self, sequence):
        '''Viterbi algorithm is a dynamic programming algorithm that allows to calculate the most probable state path for an observed sequence.
            The recurrence is given by: probability of symbol x at position t and state i
            p_i (x, t) = emission_w(x) * max {p_w(x, t-1) * p(i | w)}
        '''
        seq_len = len(sequence)    
        if seq_len == 0:
            return []
        
        viterbi = {}
        state_path = {}
        # initialize the probabilities for the first symbol
        for state in self.states:
            viterbi[state] = self.get_init_prob(state) * self.get_emission_prob(state, sequence[0])
            state_path[state] = [state]

        # compute recursively until the last element
        for t in range(1, seq_len):
            new_state_path = {}
            new_path = {}
            viterbi_tmp = {}
            for state_dest in self.states:
                intermediate_probs = []
                for state_orig in self.states:
                    prob = viterbi[state_orig] * self.get_transition_prob(state_orig, state_dest)
                    intermediate_probs.append((prob, state_orig))
                
                (max_prob, max_state) = max(intermediate_probs)
                prob = self.get_emission_prob(state_dest, sequence[t]) * max_prob      
                viterbi_tmp[state_dest] = prob
                new_state_path[state_dest] = max_state
                new_path[state_dest] = state_path[max_state] + [state_dest]
                
            viterbi = viterbi_tmp
            state_path = new_path # just keep the optimal path
    
        max_state = None
        max_prob = 0
        # among the last states find the best probability and the best path
        for state in self.states:
            if viterbi[state] > max_prob:
                max_prob = viterbi[state]
                max_state = state
                
        return (max_prob, state_path[max_state])
        
    
    def baum_welch(self, sequence):
        '''Computes an update of the emission and transition probabilities based on observed sequence
        Expectation phase: expected emission and transition probabilities are calculated based on the gamma and xi formulas
        Maximization phase: updates to the emission and transition probabilities are made based on e_hat and T_hat formulas
        '''
        print (sequence)
        seq_len = len(sequence)    
        if seq_len == 0:
            return []
        
        alpha = self.forward(sequence)
        beta = self.backward(sequence)        
        
        # Expectation phase
        # gamma: probs. for finding a symbol w at position t in state i;  product of the forward (alpha) and backward (beta) probs for all t and all i
        gamma = [{} for t in range(seq_len)] 
        # xi: probs for transition between state i at position t and state j at position t+1, for all i,j and t
        xi = [{} for t in range(seq_len - 1)]  
        for t in range(seq_len):
            # compute gamma
            sum_alpha_beta = 0
            for i in self.states:
                gamma[t][i] = alpha[t][i] * beta[t][i]
                sum_alpha_beta += gamma[t][i]
                
            for i in self.states:
                gamma[t][i] = gamma[t][i] / sum_alpha_beta
            
            # set initial probs
            if t == 0:
                self.set_init_prob(i, gamma[t][i])
            
            # compute xi values up to T - 1
            if t == seq_len - 1:
                continue
            
            sum_probs = 0
            for state_orig in self.states:
                xi[t][state_orig] = {}
                for state_dest in self.states:
                    p = alpha[t][state_orig] * self.get_transition_prob(state_orig, state_dest) * self.get_emission_prob(state_dest, sequence[t + 1]) * beta[t + 1][state_dest]
                    xi[t][state_orig][state_dest] = p 
                    sum_probs += p
                    
            for state_orig in self.states:
                for state_dest in self.states:
                    xi[t][state_orig][state_dest] /= sum_probs 
                    
        # Maximization step: with gamma and xi calculated re-estimate emissions and transitions
        # re-estimate emissions
        for i in self.states:
            denominator = 0
            for t in range(seq_len):
                denominator += gamma[t][i]

            for w in self.symbols:
                numerator = 0.0
                for t in range(seq_len):
                    if sequence[t] == w:
                        numerator += gamma[t][i]
                if denominator > 0:
                    self.set_emission_prob(i, w, numerator / denominator)
                else:
                    self.set_emission_prob(i, w, 0.0)
                
        # re-estimate transitions    
        # now that we have gamma and xi let us re-estimate
        for i in self.states:
            for j in self.states:
                denominator = 0.0
                for t in range(seq_len -1):
                    denominator += gamma[t][i]
                numerator = 0.0
                for t in range(seq_len -1):
                    numerator += xi[t][i][j]
                self.set_transition_prob(i, j, numerator / denominator)
                if denominator > 0:
                    self.set_transition_prob(i, j, numerator / denominator)
                else:
                    self.set_transition_prob(i, j, 0.0)

def test():
    initial_probs = {"5": 0.8,"M": 0.15,"3": 0.05}    
    emission_probs = {"5" : {"A": 0.20, "C": 0.30, "G":0.30, "T":0.20}, "M" : {"A":0.25, "C":0.25, "G":0.25, "T":0.25}, "3": {"A":0.35, "C":0.15, "G":0.15, "T":0.35}}
    transition_probs = {"5":{"5": 0.8,"M": 0.2,"3": 0.0},"M":{"5": 0.0,"M": 0.9,"3": 0.1}, "3":{"5": 0.0,"M": 0.0,"3": 1.0}}
    
    hmm = HiddenMarkovModel(initial_probs, emission_probs, transition_probs)
    
    # joint probability of observed sequence and state path
    sequence = "ATGCAATGCGCATGCTAAAA"
    statepath = "555555MMMMMM33333333"
    prob = hmm.joint_probability(sequence, statepath)
    print ("\nJoint probability of " + sequence + " and  " + statepath + " : " + str(prob))
    
    # forward probabilities
    seq = "ATGTGTGCACGCACCGTGCGACGCGTCGCGGAAGCTGTTATA"
    alpha = hmm.forward(seq)
    alpha_prob = sum(alpha[-1].values())
    print ("\nForward probability of " + seq + " : " + str(alpha_prob))
    
    # backward probabilities
    print ("Calculate backward probabilities")
    beta = hmm.backward(seq)

    # Optimal path with Viterbi
    seq = "ACAATGCCGTCTCCGCGACGCCTTTAATTAT"
    (probs, path) = hmm.viterbi(seq)
    print ("\nRunning Viterbi algorithm for seq " + seq)
    print ("Optimal probability: " + str(probs))
    print ("Optimal statepath: " + " ".join(path))
   

    # learn with Forward-Backward
    print ("\n\nProbabilities before learning")
    print ("Emission probabilities")
    print (hmm.emission_probs)
    print ("Transition probabilities")
    print (hmm.transition_probs)
    
    print ("Optimizing probabilities with Baum-Welch: " )
    seq = "AGGGACGCTAAGCTCGCGCGAGCGACGCCATTATAGCGTAGCTTTTTAT"
    print ("Input sequence: " + seq)
    hmm.baum_welch(seq)
    seq = "ATGTGGCGCGCGGAAGCTGTTATA"
    print ("Input sequence: " + seq)
    hmm.baum_welch(seq)
    seq = "AATCGCGAGCGGCCCGCGAAGCTGTTTTTTAATA"
    print ("Input sequence: " + seq)
    hmm.baum_welch(seq)
    seq = "ATGATGCGCTCGATGCTATCGCGCCGCGCGCGAGCGGCCCGCGAAGCTGTTTTAGTTAATAATGATATTGTA"
    print ("Input sequence: " + seq)
    hmm.baum_welch(seq)
   
    print ("Probabilities after learning")
    print ("Emission probabilities")
    print (hmm.emission_probs)
    print ("Transition probabilities")
    print (hmm.transition_probs)
    
    
test()
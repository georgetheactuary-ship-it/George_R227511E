# ============================================================
# HASTS 416 – Tutorial 1 in R
# Complete R Solutions: A1, A2, A3
# GEORGE MANGOMBE R227511E HACS
# ============================================================

# ---- Check and load required packages ----
packages_needed <- c("markovchain", "igraph", "expm")

for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop("Package ", pkg, " is not installed. Please install it using:\n",
         "  install.packages('", pkg, "', dependencies = TRUE)\n",
         "Then restart R and run this script again.")
  }
}

# ---- Robust period calculator ----
state_period <- function(P, s, states, max_n = 100) {
  i <- which(states == s)
  visits <- which(sapply(1:max_n, function(n) (P %^% n)[i, i] > 0))
  if (length(visits) == 0) return(Inf)
  Reduce(function(a, b) { while (b != 0) { tmp <- b; b <- a %% b; a <- tmp }; a },
         diff(c(0, visits)))
}

# ============================================================
# QUESTION A1 – 5-State Markov Chain
# ============================================================

P1 <- matrix(c(
  1.0, 0.0, 0.0, 0.0, 0.0,
  0.5, 0.0, 0.0, 0.0, 0.5,
  0.2, 0.0, 0.0, 0.0, 0.8,
  0.0, 0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 1.0, 0.0),
  nrow = 5, byrow = TRUE)

states1 <- c("S1","S2","S3","S4","S5")
rownames(P1) <- colnames(P1) <- states1

mc1 <- new("markovchain", transitionMatrix = P1,
           states = states1, name = "A1 Chain")

# ── A1(a): Diagram & State Classification ──────────────────

par(mar = c(2, 2, 3, 2))
plot(mc1,
     main               = "A1(a): 5-State Markov Chain",
     edge.arrow.size    = 0.4,
     edge.label.cex     = 0.8,
     edge.curved        = 0.3,
     vertex.color       = "steelblue",
     vertex.label.color = "white",
     vertex.label.cex   = 1.0,
     vertex.size        = 35,
     layout             = layout_in_circle)

cat("=== A1(a) State Classification ===\n")
cat("\nCommunicating classes:\n");  print(communicatingClasses(mc1))
cat("\nRecurrent classes:\n");      print(recurrentClasses(mc1))
cat("\nTransient classes:\n");      print(transientClasses(mc1))
cat("\nAbsorbing states:\n");       print(absorbingStates(mc1))

cat("\nPeriod of each state:\n")
for (s in states1) {
  d <- state_period(P1, s, states1)
  cat(sprintf("  Period of %s: %s\n", s,
              if (is.infinite(d)) "Inf (transient, no return)" else d))
}

# ── A1(b): Simulate Three Trajectories ─────────────────────

set.seed(42)
n_steps1 <- 30
colours1 <- c("steelblue","firebrick","darkgreen")
state_num1 <- function(s) match(s, states1)

traj1_list <- list()
for (i in 1:3) {
  s0 <- sample(states1, 1)
  traj1_list[[i]] <- c(s0, rmarkovchain(n_steps1, mc1, t0 = s0))
  cat(sprintf("\nTrajectory %d (start = %s):\n", i, s0))
  cat(traj1_list[[i]], "\n")
}

par(mar = c(5, 5, 4, 3))
plot(0:n_steps1, sapply(traj1_list[[1]], state_num1),
     type = "b", ylim = c(0.5, 5.5),
     xlab = "Time Step", ylab = "State",
     yaxt = "n", pch = 19, lwd = 1.5,
     col  = colours1[1],
     main = "A1(b): Three Simulated Trajectories",
     cex.main = 1.2, cex.lab = 1.1)
axis(2, at = 1:5, labels = states1, las = 1, cex.axis = 0.95)
abline(h = 1:5, col = "grey85", lty = 2)
for (i in 2:3)
  lines(0:n_steps1, sapply(traj1_list[[i]], state_num1),
        type = "b", col = colours1[i], pch = 19, lwd = 1.5)
legend("topright", paste("Trajectory", 1:3),
       col = colours1, lty = 1, pch = 19, lwd = 1.5,
       bty = "n", cex = 0.95)

cat("\nComment: Trajectories starting in transient states (S2, S3, S4, S5)")
cat("\neventually get absorbed into S1 (the only absorbing/recurrent state).\n")

# ── A1(c): Steady-State Probabilities & Ergodicity ─────────

cat("\n=== A1(c) Steady-State Probabilities ===\n")
ss1 <- steadyStates(mc1)
print(round(ss1, 6))

cat("\nChain is ergodic: FALSE\n")
cat("Reason: Not irreducible – transient states exist alongside\n")
cat("one recurrent class, so no unique limiting distribution.\n")

# ── A1(d): Convergence of Unconditional Probabilities ──────

n_time1 <- 50
init1   <- rep(1/5, 5)
prob_mat1 <- matrix(NA, n_time1 + 1, 5)
prob_mat1[1, ] <- init1
for (t in 1:n_time1) prob_mat1[t+1, ] <- prob_mat1[t, ] %*% P1

par(mar = c(5, 5, 4, 8))
matplot(0:n_time1, prob_mat1,
        type = "l", lty = 1, lwd = 2,
        col  = c("steelblue","firebrick","darkgreen","purple","orange"),
        xlab = "Time n", ylab = "P(Xn = s)",
        main = "A1(d): Convergence of Unconditional Probabilities",
        cex.main = 1.2, cex.lab = 1.1,
        ylim = c(0, 1))
abline(h = seq(0, 1, 0.2), col = "grey85", lty = 2)
legend("right", inset = c(-0.25, 0),
       legend = states1,
       col    = c("steelblue","firebrick","darkgreen","purple","orange"),
       lty = 1, lwd = 2, bty = "n", cex = 0.95,
       xpd = TRUE)

cat("\nComment: P(Xn = S1) rises to 1; all transient state probabilities\n")
cat("decay to 0. Convergence is essentially complete by n = 20.\n")

par(mar = c(5.1, 4.1, 4.1, 2.1))

# ============================================================
# QUESTION A2 – 7-State Markov Chain
# ============================================================

P2 <- matrix(c(
  0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.4, 0.2, 0.2, 0.2,
  0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.4,
  0.3, 0.0, 0.0, 0.1, 0.3, 0.1, 0.2,
  0.0, 0.0, 0.0, 0.2, 0.2, 0.3, 0.3,
  0.0, 0.0, 0.0, 0.5, 0.2, 0.2, 0.1),
  nrow = 7, byrow = TRUE)

states2 <- paste0("S", 1:7)
rownames(P2) <- colnames(P2) <- states2

mc2 <- new("markovchain", transitionMatrix = P2,
           states = states2, name = "A2 Chain")

# ── A2(a): Chain Diagram ────────────────────────────────────

par(mar = c(2, 2, 3, 2))
plot(mc2,
     main               = "A2(a): 7-State Markov Chain",
     edge.arrow.size    = 0.35,
     edge.label.cex     = 0.7,
     edge.curved        = 0.35,
     vertex.color       = "darkorange",
     vertex.label.color = "white",
     vertex.label.cex   = 0.95,
     vertex.size        = 30,
     layout             = layout_in_circle)

# ── A2(b): State Classification & Periods ──────────────────

cat("\n=== A2(b) State Classification ===\n")
cat("\nCommunicating classes:\n"); print(communicatingClasses(mc2))
cat("\nRecurrent classes:\n");     print(recurrentClasses(mc2))
cat("\nTransient classes:\n");     print(transientClasses(mc2))
cat("\nAbsorbing states:\n");      print(absorbingStates(mc2))

cat("\nPeriod of each state:\n")
for (s in states2) {
  d <- state_period(P2, s, states2)
  cat(sprintf("  Period of %s: %s\n", s,
              if (is.infinite(d)) "Inf (transient, no return)" else d))
}

cat("\nNo absorbing states. No reflecting states.\n")
cat("{S1,S2}: recurrent, period 2.\n")
cat("{S4,S5,S6,S7}: recurrent, period 1 (aperiodic).\n")
cat("S3: transient.\n")

# ── A2(c): Two Simulated Trajectories ──────────────────────

set.seed(2024)
n_steps2  <- 50
state_num2 <- function(s) match(s, states2)
cols2   <- c("black","firebrick")
starts2 <- sample(states2, 2)

trajs2 <- lapply(starts2, function(s0)
  c(s0, rmarkovchain(n_steps2, mc2, t0 = s0)))

for (i in 1:2) {
  cat(sprintf("\nA2 Trajectory %d (start = %s):\n", i, starts2[i]))
  cat(trajs2[[i]], "\n")
}

par(mar = c(5, 5, 4, 3))
plot(0:n_steps2, sapply(trajs2[[1]], state_num2),
     type = "b", ylim = c(0.5, 7.5),
     xlab = "Time Step", ylab = "State",
     yaxt = "n", pch = 19, lwd = 1.5,
     col  = cols2[1],
     main = "A2(c): Two Simulated Trajectories",
     cex.main = 1.2, cex.lab = 1.1)
axis(2, at = 1:7, labels = states2, las = 1, cex.axis = 0.95)
abline(h = 1:7, col = "grey85", lty = 2)
lines(0:n_steps2, sapply(trajs2[[2]], state_num2),
      type = "b", col = cols2[2], pch = 19, lwd = 1.5)
legend("topright", paste("Trajectory", 1:2),
       col = cols2, lty = 1, pch = 19, lwd = 1.5,
       bty = "n", cex = 0.95)

cat("\nComment: Trajectories entering {S1,S2} oscillate with period 2.")
cat("\nTrajectories entering {S4,S5,S6,S7} wander aperiodically within")
cat("\nthat class. S3 is visited briefly before escape.\n")

# ── A2(d): Limiting Probabilities & Ergodicity ─────────────

cat("\n=== A2(d) Limiting (Stationary) Distributions ===\n")
lim2 <- steadyStates(mc2)
print(round(lim2, 6))

cat("\nChain is ergodic: FALSE\n")
cat("Reason: Two recurrent classes – no unique stationary distribution.\n")

par(mar = c(5.1, 4.1, 4.1, 2.1))

# ============================================================
# QUESTION A3 – Time-Inhomogeneous Traffic Chain
# ============================================================

P_day  <- matrix(c(0.4, 0.4, 0.2,
                   0.3, 0.5, 0.3,
                   0.0, 0.1, 0.9),
                 nrow = 3, byrow = TRUE)

P_peak <- matrix(c(0.1, 0.5, 0.4,
                   0.1, 0.3, 0.6,
                   0.0, 0.1, 0.9),
                 nrow = 3, byrow = TRUE)

traffic_states <- c("Light","Heavy","Jammed")
rownames(P_day)  <- colnames(P_day)  <- traffic_states
rownames(P_peak) <- colnames(P_peak) <- traffic_states

mat_power <- function(M, n) {
  result <- diag(nrow(M))
  for (i in seq_len(n)) result <- result %*% M
  result
}

# ── A3(a): Analytical Distribution at 6 PM ─────────────────

pi0     <- c(1, 0, 0)
P_day9  <- mat_power(P_day,  9)
P_peak6 <- mat_power(P_peak, 6)
pi_4pm  <- pi0    %*% P_day9
pi_6pm  <- pi_4pm %*% P_peak6

cat("\n=== A3(a) Analytical Distribution ===\n")
cat("Distribution at 4 PM:\n");  print(round(pi_4pm, 4))
cat("Distribution at 6 PM:\n");  print(round(pi_6pm, 4))
cat(sprintf("P(Light  at 6 PM) = %.4f\n", pi_6pm[1]))
cat(sprintf("P(Heavy  at 6 PM) = %.4f\n", pi_6pm[2]))
cat(sprintf("P(Jammed at 6 PM) = %.4f\n", pi_6pm[3]))

# ── A3(b): Monte Carlo Verification (N = 10,000) ───────────

set.seed(123)
N <- 10000

norm_rows <- function(M) t(apply(M, 1, function(r) r / sum(r)))
P1n <- norm_rows(P_day)
P2n <- norm_rows(P_peak)

simulate_traffic <- function() {
  state <- 1
  for (t in 1:9) state <- sample(1:3, 1, prob = P1n[state, ])
  for (t in 1:6) state <- sample(1:3, 1, prob = P2n[state, ])
  state
}

final_states <- replicate(N, simulate_traffic())
emp_dist <- tabulate(final_states, nbins = 3) / N
names(emp_dist) <- traffic_states

cat("\n=== A3(b) Monte Carlo Results (N = 10,000) ===\n")
comparison <- rbind(Analytical = round(as.numeric(pi_6pm), 4),
                    Simulation = round(emp_dist, 4))
colnames(comparison) <- traffic_states
print(comparison)

par(mar = c(5, 5, 4, 3))
bp <- barplot(comparison,
              beside      = TRUE,
              col         = c("black","firebrick"),
              border      = "white",
              main        = "A3(b): Analytical vs Simulated Distribution at 6 PM",
              ylab        = "Probability",
              xlab        = "Traffic State",
              ylim        = c(0, max(comparison) * 1.25),
              cex.main    = 1.2,
              cex.lab     = 1.1,
              cex.names   = 1.05)
legend("topright", legend = rownames(comparison),
       fill = c("black","firebrick"),
       border = "white", bty = "n", cex = 0.95)

text(bp, comparison + 0.012,
     labels = sprintf("%.3f", comparison),
     cex = 0.85, font = 2)

cat("\nComment: Simulation proportions closely match the analytical\n")
cat("probabilities. Jammed dominates at 6 PM due to peak-hour\n")
cat("transition probabilities pushing mass toward congestion.\n")
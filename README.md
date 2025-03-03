# Track Heatmap Generator

## 📌 Overview
This project provides a **Track Heatmap Generator** that visualizes GPS track data in a **grid-based format**. The program processes **GPS trackpoints** from a series of recorded tracks and generates a **heatmap**, representing the density of trackpoints within geographic grid cells.

This heatmap helps analyze **frequently traveled areas**, **popular paths**, and **coverage density** over a specific region.

---

## 🚀 Features
- **Processes GPS Track Data** from multiple recorded segments.
- **Computes a grid-based heatmap** based on trackpoint density.
- **Supports configurable cell width and height** for different levels of granularity.
- **Handles arbitrary GPS locations and computes coverage efficiently**.
- **Efficient allocation & memory management**, ensuring scalability.

---

## 🛠 Installation & Setup
### **Requirements**
- **C Compiler (gcc, clang, etc.)**
- **Make (optional, for easy compilation)**
- **Math Library (`-lm` needed for compilation)**

### **Compiling the Program**
To compile the program, run:
```bash
gcc -o heatmap track.c trackpoint.c location.c -lm

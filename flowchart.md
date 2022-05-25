graph TD
    A(START) --> B[/"1️⃣ &nbsp;&nbsp;&nbsp Load <b>points</b><br> and <b>target</b> variable"/]
    B --> C["2️⃣ &nbsp;&nbsp; Compute <b>distance</b> matrix"]
    C --> D["3️⃣ &nbsp;&nbsp; Compute <b>variability</b> matrix"]
    D --> G{"4️⃣ &nbsp;&nbsp; Are all links <br> <b>connected</b>?"}
    G -- YES --> L["9️⃣ &nbsp;&nbsp; Compute <b>SGI</b> applying <b>Lorenz curve</b> </br> approach to variability for each h-distance"]
    G --NO--> H["5️⃣ &nbsp;&nbsp; Build <b>h-distance</b> applying <br> <b>MaxMin</b> method to <br> distance matrix"]
    H --> I["6️⃣ &nbsp;&nbsp; Build <b>contiguity</b> and <br><b>weight</b> matrix considering <br> a <b>range distance band</b> <br> between a min-threshold <br> (previous h-distance) <br> and a max-threshold (h-distance)"]
    I --> J["7️⃣ &nbsp;&nbsp; Compute <b>total variability</b> </br>multiplying the weight matrix with <br> variability matrix and summing the elements"]
    J --> K["8️⃣ &nbsp;&nbsp; Remove from distance matrix </br> the links connected in the current iteration"]
    K --> G
    L ----> M(END)
    
    

# DrugDiscov - AI-Powered Drug Discovery Platform

DrugDiscov is a web application for AI-driven drug discovery and molecular optimization. This platform enables researchers and pharmaceutical teams to create projects, generate and optimize molecules, and view detailed AI-generated reports on drug candidates.

## Technology Stack

- **Frontend**: React.js with TypeScript
- **UI Libraries**: Tailwind CSS, shadcn/ui components
- **State Management**: React Query
- **Routing**: React Router
- **Icons**: Lucide React
- **Charts & Visualizations**: Recharts (for analytics)

## Application Structure

- **Dashboard**: Overview of current projects and key metrics
- **Projects**: Management of drug discovery projects
- **Molecule Generator**: AI-based tool for generating new molecules
- **Molecule Optimization**: Fine-tuning molecules for specific properties
- **Molecule Preview**: 3D visualization of molecules and key characteristics
- **Project Reports**: Detailed AI analysis of drug candidates

## Backend Integration Requirements

### API Endpoints Needed

1. **Projects API**
   - `GET /api/projects` - List all projects
   - `GET /api/projects/:id` - Get project details
   - `POST /api/projects` - Create new project
   - `PUT /api/projects/:id` - Update project
   - `DELETE /api/projects/:id` - Delete project

2. **Molecules API**
   - `GET /api/projects/:id/molecules` - List molecules for a project
   - `GET /api/molecules/:id` - Get molecule details
   - `POST /api/molecules/generate` - Generate new molecules
   - `POST /api/molecules/optimize` - Optimize existing molecules
   - `GET /api/molecules/:id/preview` - Get 3D visualization data

3. **Reports API**
   - `GET /api/projects/:id/report` - Get project report
   - `GET /api/molecules/:id/report` - Get molecule report

### Data Structures

#### Project
```typescript
interface Project {
  id: string;
  name: string;
  description: string;
  targetProtein: string;
  status: 'Planning' | 'In Progress' | 'On Hold' | 'Completed';
  therapeuticArea: string;
  collaborators: string[];
  desiredProperties: {
    solubility: 'High' | 'Medium' | 'Low';
    bindingAffinity: 'High' | 'Medium' | 'Low';
    toxicity: 'Low' | 'Medium' | 'High';
    bioavailability: 'High' | 'Medium' | 'Low';
  };
  molecules: number;
  lastUpdated: string;
}
```

#### Molecule
```typescript
interface Molecule {
  id: string;
  projectId: string;
  name: string;
  smiles: string;
  inchi: string;
  molecularWeight: number;
  logP: number;
  hBondDonors: number;
  hBondAcceptors: number;
  rotatableBonds: number;
  properties: {
    solubility: string;
    bindingAffinity: string;
    toxicity: string;
    bioavailability: string;
    syntheticAccessibility: number;
  };
  targets: Array<{
    name: string;
    probability: number;
    confidence: 'High' | 'Medium' | 'Low';
  }>;
  created: string;
}
```

#### Report
```typescript
interface Report {
  id: string;
  subjectId: string; // Project ID or Molecule ID
  subjectType: 'project' | 'molecule';
  executiveSummary: {
    drugLikenessScore: number;
    lipinskiViolations: number;
    toxicityConcerns: number;
    syntheticAccessibility: number;
    bbbPermeability: string;
    overallClassification: string;
  };
  physicochemicalProperties: {
    // Properties as shown in the report
  };
  admetPredictions: {
    // ADMET details
  };
  toxicityProfile: {
    // Toxicity information
  };
  targetPrediction: Array<{
    // Target information
  }>;
  drugInteractions: Array<{
    // Interaction data
  }>;
  // Additional report sections
  generated: string;
}
```

### Authentication

The frontend is prepared to work with JWT-based authentication. Backend endpoints should:

1. Implement token-based authentication
2. Return appropriate headers for CORS
3. Include role-based access control for different user types

### 3D Molecule Visualization

For molecule visualization, the frontend expects:
- 3D coordinates in PDB or similar format
- Or the ability to generate 3D structure from SMILES notation

## Development Setup

1. Clone the repository
2. Install dependencies: `npm install`
3. Start the development server: `npm run dev`
4. Access the application at: `http://localhost:8080`

## Next Steps for Backend Integration

1. Implement the API endpoints described above
2. Ensure proper error handling and response formats
3. Document any changes or additional requirements
4. Set up authentication and authorization
5. Implement AI models for molecule generation and report creation

## Contact

For questions or further clarification, please reach out to the frontend development team.

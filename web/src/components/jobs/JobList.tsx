"use client";

import { useEffect, useState } from "react";
import { JobProgress, Job } from "./JobProgress";
import { RefreshCw } from "lucide-react";

const API_URL = "http://localhost:8000/api";

export function JobList() {
  const [jobs, setJobs] = useState<Job[]>([]);
  const [loading, setLoading] = useState(true);

  const fetchJobs = async () => {
    try {
      const res = await fetch(`${API_URL}/jobs`);
      if (res.ok) {
        setJobs(await res.json());
      }
    } catch {
      // API not available
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    fetchJobs();
    const interval = setInterval(fetchJobs, 5000);
    return () => clearInterval(interval);
  }, []);

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between">
        <h2 className="text-lg font-semibold">Jobs</h2>
        <button
          onClick={fetchJobs}
          className="rounded-lg p-2 text-slate-400 hover:bg-slate-800 hover:text-slate-200 transition-colors"
        >
          <RefreshCw className="h-4 w-4" />
        </button>
      </div>

      {loading ? (
        <p className="text-slate-500 text-sm">Loading jobs...</p>
      ) : jobs.length === 0 ? (
        <p className="text-slate-500 text-sm">No jobs yet. Long-running tools will appear here.</p>
      ) : (
        <div className="space-y-3">
          {jobs.map((job) => (
            <JobProgress key={job.id} job={job} />
          ))}
        </div>
      )}
    </div>
  );
}
